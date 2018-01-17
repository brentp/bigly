package main

import (
	"compress/gzip"
	"encoding/json"
	"fmt"
	"html/template"
	"io"
	"log"
	"math"
	"net/http"
	"os"
	"sort"
	"strconv"
	"strings"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/brentp/bigly"
	"github.com/brentp/faidx"
	chartjs "github.com/brentp/go-chartjs"
	"github.com/brentp/goleft/covstats"
	"github.com/brentp/xopen"
)

var MaxRegion int = 1e7

// MinSoftClips determines how many softclips are required Atoi
// each position in order to be plotted. Higher values remove
// noise.
// TODO: make this a cliarg parameter.
const MinSoftClips = 5

// MinSoftClipProportion determines the proporition of reads that must be soft clipped in order to by shown.
const MinSoftClipProportion = 0.1

type cliarg struct {
	bigly.Options
	Reference string       `arg:"-r,help:optional path to reference fasta."`
	BamPath   []string     `arg:"positional,required"`
	ref       *faidx.Faidx `arg:"-"`
	Port      int          `arg:"-p,help:server port"`
	// maps from sample id to bam path.
	paths map[string]string `arg:"-"`

	bamStats *covstats.Stats `arg:"-"`
}

func mean(arr []int32) int {
	if len(arr) == 0 {
		return 0
	}
	var s float64
	for _, a := range arr {
		s += float64(a)
	}
	return int(0.5 + s/float64(len(arr)))
}

// satisfy the required interface with this struct and methods.
type xy struct {
	x []float64
	y []float64
	r []float64
}

func (v xy) Xs() []float64 {
	return v.x
}
func (v xy) Ys() []float64 {
	return v.y
}
func (v xy) Rs() []float64 {
	return v.r
}

func (c cliarg) Version() string {
	return "bigly 0.2.0"
}

func getShortName(b string) string {

	fh, err := os.Open(b)
	if err != nil {
		log.Fatal(err)
	}
	defer fh.Close()
	br, err := bam.NewReader(fh, 2)
	if err != nil {
		log.Fatal(err)
	}
	defer br.Close()
	m := make(map[string]bool)
	for _, rg := range br.Header().RGs() {
		m[rg.Get(sam.Tag([2]byte{'S', 'M'}))] = true
	}
	if len(m) > 1 {
		log.Println("warning: more than one tag for %s", b)
	}
	for sm := range m {
		return sm
	}
	return ""
}

func parseRegion(region string) (chrom string, start, end int, err error) {
	chromse := strings.Split(region, ":")
	if len(chromse) != 2 {
		return "", 0, 0, fmt.Errorf("expected a region like {chrom}:{start}-{end}. Got %s", region)
	}
	se := strings.Split(chromse[1], "-")
	start, err = strconv.Atoi(strings.Replace(se[0], ",", "", -1))
	if err != nil {
		return "", 0, 0, fmt.Errorf("unable to parse region %s", region)
	}
	end, err = strconv.Atoi(strings.Replace(se[1], ",", "", -1))
	if err != nil {
		return "", 0, 0, fmt.Errorf("unable to parse region %s", region)
	}
	return chromse[0], start - 1, end, nil
}

type ifill struct {
	SampleNames []string
	Region      string
}

func getRegion(region string) string {
	if region == "" {
		region = "7:71211015-71213751"
	}
	return region
}

func (cli *cliarg) ServeIndex(w http.ResponseWriter, r *http.Request) {
	t, err := template.ParseFiles("templates/chartjs.tmpl")
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err)
		log.Fatal(err)
	}
	cli.paths = make(map[string]string, len(cli.BamPath))
	samples := make([]string, 0, len(cli.BamPath))
	for _, p := range cli.BamPath {
		name := getShortName(p)
		cli.paths[name] = p
		samples = append(samples, name)
	}
	if err = t.Execute(w, ifill{SampleNames: samples, Region: getRegion(r.FormValue("region"))}); err != nil {
		log.Fatal(err)
	}
	wtr, _ := xopen.Wopen("index.html")
	t.Execute(wtr, cli)
	wtr.Close()
}

type tfill struct {
	Depths    xy
	Splitters xy
	Inserts   xy
	Softs     xy
}

func abs(p float64) float64 {
	if p < 0 {
		return -p
	}
	return p
}

func toPct(a, depth float64) float64 {
	return a
	if depth == 0 {
		if a > 0 {
			return 100.0
		}
		return 0
	}
	return 100.0 * float64(a) / float64(depth)
}

func fmax(a, b uint32) float64 {
	v := b
	if a > b {
		v = a
	}
	if v > 1000 {
		v = 1000
	}
	return float64(v)
}

// sometimes we get a single point with a signal for a large insert size when
// there is an island of coverage in a no coverage area.
// this removes large (2e4) values of insert size that are preceded by no coverage.
func removeSoleOutliers(inserts xy) {
	for i := 1; i < len(inserts.y); i++ {
		if inserts.y[i-1] == 0 && inserts.y[i]-inserts.y[i-1] > 2e4 && i < len(inserts.y)-1 && inserts.y[i+1] == 0 {
			inserts.y[i] = 0
		}
	}
}

func (cli *cliarg) ServeHTTP(w http.ResponseWriter, r *http.Request) {
	ctx := r.Context()
	w.Header().Set("Content-Type", "application/json")
	w.Header().Set("Content-Encoding", "gzip")
	po := strings.Split(strings.TrimSpace(r.URL.Path[len("/data/"):]), "/")
	if len(po) != 2 {
		http.Error(w, fmt.Sprintf("expected a path like {bam}/{chrom}:{start}-{end}. Got %s", r.URL.Path[1:]), http.StatusInternalServerError)
		fmt.Fprintln(os.Stderr, fmt.Sprintf("expected a path like {bam}/{chrom}:{start}-{end}. Got %s", r.URL.Path[1:]))
		fmt.Fprintln(os.Stderr, fmt.Sprintf("Got %v\n", po))
		return
	}
	name, region := po[0], po[1]
	log.Println("got request:", name, region)

	chrom, start, end, err := parseRegion(region)
	if end-start > MaxRegion {
		return
	}
	if err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		fmt.Fprintln(os.Stderr, err.Error())
		return
	}

	bamPath := cli.paths[name]

	it := bigly.Up(bamPath, cli.Options, bigly.Position{Chrom: chrom, Start: start, End: end}, cli.ref)
	tf := tfill{Depths: xy{}, Splitters: xy{}, Inserts: xy{}, Softs: xy{}}
	tf.Inserts.x = append(tf.Inserts.x, float64(start))
	tf.Inserts.y = append(tf.Inserts.y, math.NaN())

	splits := make(map[int]int)
	for it.Next() {
		p := it.Pile()
		select {
		case <-ctx.Done():
			// got cancelled by the client.
			return
		default:
		}

		if len(tf.Depths.y) == 0 || abs(tf.Depths.y[len(tf.Depths.y)-1]-float64(p.Depth)) >= 1 {
			if len(tf.Depths.y) != 0 && tf.Depths.y[len(tf.Depths.y)-1] == 0 {
				tf.Depths.y = append(tf.Depths.y, 0)
				tf.Depths.x = append(tf.Depths.x, float64(p.Pos-1))
			}
			tf.Depths.y = append(tf.Depths.y, float64(p.Depth))
			tf.Depths.x = append(tf.Depths.x, float64(p.Pos))
		}

		if p.SoftStarts+p.SoftEnds >= MinSoftClips && float64(p.SoftStarts+p.SoftEnds)/float64(p.Depth) > MinSoftClipProportion {
			tf.Softs.x = append(tf.Softs.x, float64(p.Pos))
			tf.Softs.x = append(tf.Softs.x, float64(p.Pos))
			tf.Softs.x = append(tf.Softs.x, float64(p.Pos))

			tf.Softs.y = append(tf.Softs.y, 0)
			tf.Softs.y = append(tf.Softs.y, float64(p.SoftStarts+p.SoftEnds))
			tf.Softs.y = append(tf.Softs.y, math.NaN())
		}

		if p.Splitters > 1 {
			posns := getPositions(p.SplitterPositions)
			m, c := bigly.Mode(posns)
			if c > 1 && m+100 > start && m-100 < end {
				if _, ok := splits[m]; !ok {
					splits[m] = 0
				}
				splits[m] += c
			}
		}
		// we take the max of left and right as it gives a cleaner signal than mean.
		in := fmax(uint32(mean(p.InsertSizeLPs)), uint32(mean(p.InsertSizeRMs)))
		last := tf.Inserts.y[len(tf.Inserts.y)-1]
		if last == 0 || math.IsNaN(last) {
			tf.Inserts.y = append(tf.Inserts.y, math.NaN())
			tf.Inserts.y = append(tf.Inserts.y, math.NaN())

			tf.Inserts.x = append(tf.Inserts.x, tf.Inserts.x[len(tf.Inserts.x)-1])
			tf.Inserts.x = append(tf.Inserts.x, float64(p.Pos-1))
		}
		tf.Inserts.x = append(tf.Inserts.x, float64(p.Pos))
		tf.Inserts.y = append(tf.Inserts.y, float64(in))
	}
	//removeSoleOutliers(tf.Inserts)
	max := 1
	sites := make([]int, 0, len(splits))
	for pos, v := range splits {
		sites = append(sites, pos)
		if v > max {
			max = v
		}
	}
	sort.Ints(sites)

	for k, pos := range sites {
		v := splits[pos]
		if k == 0 {
			tf.Splitters.x = append(tf.Splitters.x, float64(pos)-1)
			tf.Splitters.y = append(tf.Splitters.y, 0)
		} else if float64(pos)-tf.Splitters.x[k-1] > 1 {
			tf.Splitters.x = append(tf.Splitters.x, tf.Splitters.x[len(tf.Splitters.x)-1])
			tf.Splitters.y = append(tf.Splitters.y, 0)
			tf.Splitters.x = append(tf.Splitters.x, float64(pos)-1)
			tf.Splitters.y = append(tf.Splitters.y, math.NaN())
		}
		tf.Splitters.x = append(tf.Splitters.x, float64(pos))
		tf.Splitters.y = append(tf.Splitters.y, toPct(float64(v), float64(max)))
	}
	if len(sites) > 1 {
		tf.Splitters.x = append(tf.Splitters.x, tf.Splitters.x[len(tf.Splitters.x)-1]+1)
		tf.Splitters.y = append(tf.Splitters.y, math.NaN())
	}

	if err := it.Error(); err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		log.Println(err)
		return
	}
	it.Close()
	if err := writeChart(w, tf, start, end); err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		log.Println(err)
	}
}

func getPositions(ps []bigly.Position) []int {
	m := make([]int, 0, 2*len(ps))
	for _, p := range ps {
		m = append(m, p.Start, p.End)
	}
	return m
}

func writeChart(w io.Writer, tf tfill, start, end int) error {
	chart := chartjs.Chart{Label: "bigly-chart"}
	xtick := &chartjs.Tick{Min: float64(start), Max: float64(end)}
	right2, err := chart.AddYAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Right})
	if err != nil {
		return err
	}
	right1, err := chart.AddYAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Right})
	if err != nil {
		return err
	}
	left2, err := chart.AddYAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Left})
	if err != nil {
		return err
	}
	left1, err := chart.AddYAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Left})
	if err != nil {
		return err
	}

	if _, err = chart.AddXAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Bottom, Display: chartjs.True, ScaleLabel: &chartjs.ScaleLabel{LabelString: "genomic position", Display: chartjs.True}, Tick: xtick}); err != nil {
		return err
	}

	chart.AddDataset(chartjs.Dataset{
		Data: tf.Splitters, Label: "splitters", Type: chartjs.Line, YAxisID: right2,
	})

	chart.AddDataset(chartjs.Dataset{
		Data: tf.Softs, Label: "soft-clips", Type: chartjs.Line, YAxisID: right1,
	})

	chart.AddDataset(chartjs.Dataset{
		Data: tf.Inserts, Label: "insert-size", Type: chartjs.Line, YAxisID: left2,
	})

	chart.AddDataset(chartjs.Dataset{
		Data: tf.Depths, Label: "depth", Type: chartjs.Line, YAxisID: left1,
	})
	chart.Options.Responsive = chartjs.True
	chart.Options.MaintainAspectRatio = chartjs.False

	//buf, err := json.MarshalIndent(chart, "\n", "  ")
	buf, err := json.Marshal(chart)
	if err != nil {
		return err
	}

	json := fmt.Sprintf(`{"chart": %s}`, string(buf))

	gz := gzip.NewWriter(w)
	defer gz.Close()
	_, err = gz.Write([]byte(json))
	if err != nil {
		return err
	}

	return nil
}

func main() {

	cli := &cliarg{}
	cli.Options.MinBaseQuality = 10
	cli.Options.ConcordantCutoff = 10000
	cli.Options.MinMappingQuality = 5
	cli.Options.MinClipLength = 15
	cli.Options.SplitterVerbosity = 1
	cli.Port = 5000
	arg.MustParse(cli)
	if cli.ExcludeFlag == 0 {
		cli.ExcludeFlag = uint16(sam.Unmapped | sam.QCFail | sam.Duplicate)
	}

	// a path like /data/sample/1:1234-5678
	http.HandleFunc("/data/", cli.ServeHTTP)

	// fills the template
	http.HandleFunc("/", cli.ServeIndex)

	if e := http.ListenAndServe(fmt.Sprintf(":%d", cli.Port), nil); e != nil {
		log.Fatal(e)
	}

}
