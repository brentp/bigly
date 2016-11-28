package main

import (
	"compress/gzip"
	"encoding/json"
	"fmt"
	"html/template"
	"io"
	"log"
	"net/http"
	"os"
	"strconv"
	"strings"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/hts/sam"
	"github.com/brentp/bigly"
	"github.com/brentp/faidx"
	chartjs "github.com/brentp/go-chartjs"
	"github.com/brentp/vcfgo"
	"github.com/brentp/xopen"
)

type cliarg struct {
	bigly.Options
	Reference string       `arg:"-r,help:optional path to reference fasta."`
	VCF       string       `arg:"-v,help:optional vcf with variants to examine."`
	BamPath   []string     `arg:"positional,required"`
	ref       *faidx.Faidx `arg:"-"`
	Port      int          `arg:"-p,help:server port"`
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
	rgs := strings.Split(b, "/")
	rg := rgs[len(rgs)-1]
	return rg
}

func parseRegion(region string) (chrom string, start, end int, err error) {
	chromse := strings.Split(region, ":")
	if len(chromse) != 2 {
		return "", 0, 0, fmt.Errorf("expected a region like {chrom}:{start}-{end}. Got %s", region)
	}
	se := strings.Split(chromse[1], "-")
	start, err = strconv.Atoi(se[0])
	if err != nil {
		return "", 0, 0, fmt.Errorf("unable to parse region %s", region)
	}
	end, err = strconv.Atoi(se[1])
	if err != nil {
		return "", 0, 0, fmt.Errorf("unable to parse region %s", region)
	}
	return chromse[0], start - 1, end, nil
}

type ifill struct {
	BamPath []string
	Regions []string
	// contains, e.. 0/1, 0/0, from the VCF for the genotype of each sample.
	Genotypes [][]string
}

func getRegions(cli cliarg) ([]string, [][]string) {
	if cli.VCF == "" {
		return []string{}, nil
	}
	genotypes := make([][]string, 0, 20)
	m := make([]string, 0, 20)
	rdr, err := xopen.Ropen(cli.VCF)
	if err != nil {
		panic(err)
	}
	vcf, err := vcfgo.NewReader(rdr, false)
	defer rdr.Close()
	if err != nil {
		panic(err)
	}
	defer vcf.Close()
	for {
		rec := vcf.Read()
		if rec == nil {
			break
		}
		ex := uint32(float64(rec.End()-rec.Start()) / 5.0)
		m = append(m, fmt.Sprintf("%s:%d-%d", rec.Chrom(), rec.Start()-ex, rec.End()+ex))

		gts := make([]string, 0, len(rec.Samples))
		for _, s := range rec.Samples {
			gts = append(gts, s.Fields["GT"])
		}
		genotypes = append(genotypes, gts)

	}
	if e := vcf.Error(); e != nil {
		fmt.Fprintln(os.Stderr, e.Error())
	}
	return m, genotypes
}

// TODO: how to map between sample names in VCF and bam files.
func (cli cliarg) ServeIndex(w http.ResponseWriter, r *http.Request) {
	t, err := template.ParseFiles("templates/chartjs.tmpl")
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err)
		log.Fatal(err)
	}
	paths := make([]string, len(cli.BamPath))
	for i, p := range cli.BamPath {
		paths[i] = getShortName(p)
	}
	regs, gts := getRegions(cli)
	if err = t.Execute(w, ifill{BamPath: paths, Regions: regs, Genotypes: gts}); err != nil {
		log.Fatal(err)
	}
	wtr, _ := xopen.Wopen("index.html")
	t.Execute(wtr, cli)
	wtr.Close()
}

func getBamPath(cli cliarg, name string) string {
	for _, p := range cli.BamPath {
		if getShortName(p) == name {
			return p
		}
	}
	return ""
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

func mustMarshal(v []xy) string {
	s, e := json.Marshal(v)
	if e != nil {
		panic(e)
	}
	return string(s)
}

func (cli cliarg) ServeHTTP(w http.ResponseWriter, r *http.Request) {
	w.Header().Set("Content-Type", "application/json")
	w.Header().Set("Content-Encoding", "gzip")
	po := strings.Split(strings.TrimSpace(r.URL.Path[len("/data/"):]), "/")
	if len(po) != 2 {
		http.Error(w, fmt.Sprintf("expected a path like {bam}/{chrom}:{start}-{end}. Got %s", r.URL.Path[1:]), http.StatusInternalServerError)
		fmt.Fprintln(os.Stderr, fmt.Sprintf("expected a path like {bam}/{chrom}:{start}-{end}. Got %s", r.URL.Path[1:]))
		return
	}
	name, region := po[0], po[1]

	chrom, start, end, err := parseRegion(region)
	if err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		fmt.Fprintln(os.Stderr, err.Error())
		return
	}

	bamPath := getBamPath(cli, name)

	it := bigly.Up(bamPath, cli.Options, bigly.Position{chrom, start, end}, cli.ref)
	tf := tfill{xy{}, xy{}, xy{}, xy{}}
	splits := make(map[int]int)
	for it.Next() {
		p := it.Pile()
		if len(tf.Depths.y) == 0 || abs(tf.Depths.y[len(tf.Depths.y)-1]-float64(p.Depth)) >= 1 {
			if len(tf.Depths.y) != 0 && tf.Depths.y[len(tf.Depths.y)-1] == 0 {
				tf.Depths.y = append(tf.Depths.y, 0)
				tf.Depths.x = append(tf.Depths.x, float64(p.Pos-1))
			}
			tf.Depths.y = append(tf.Depths.y, float64(p.Depth))
			tf.Depths.x = append(tf.Depths.x, float64(p.Pos))
		}
		if p.SoftStarts > 1 || p.SoftEnds > 1 {
			tf.Softs.x = append(tf.Softs.x, float64(p.Pos))
			tf.Softs.y = append(tf.Softs.y, float64(p.SoftStarts+p.SoftEnds))
		}
		if p.Splitters > 1 {
			m, c := bigly.Mode(p.SplitterPositions)
			if c > 1 {
				if _, ok := splits[m]; !ok {
					splits[m] = 0
				}
				splits[m] += c
			}
		}
		in := float64(p.MeanInsertSizeLP + p.MeanInsertSizeRM)
		if len(tf.Inserts.x) == 0 || (len(tf.Inserts.x) > 0 && in != tf.Inserts.y[len(tf.Inserts.y)-1]) {
			if len(tf.Inserts.y) > 0 && tf.Inserts.y[len(tf.Inserts.y)-1] == 0 {
				tf.Inserts.y = append(tf.Inserts.y, 0)
				tf.Inserts.x = append(tf.Inserts.x, float64(p.Pos-1))
			}
			tf.Inserts.x = append(tf.Inserts.x, float64(p.Pos))
			tf.Inserts.y = append(tf.Inserts.y, float64(in))
		}
	}
	max := 1
	for _, v := range splits {
		if v > max {
			max = v
		}
	}

	for k, v := range splits {
		tf.Splitters.x = append(tf.Splitters.x, float64(k))
		tf.Splitters.y = append(tf.Splitters.y, toPct(float64(v), float64(max)))
	}

	if err := it.Error(); err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		log.Println(err)
		return
	}
	it.Close()
	if err := writeChart(w, tf); err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		log.Println(err)
	}
}

func writeChart(w io.Writer, tf tfill) error {
	chart := chartjs.Chart{Label: "bigly-chart"}
	right2, err := chart.AddYAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Right,
		ScaleLabel: &chartjs.ScaleLabel{LabelString: "splitters", Display: chartjs.True}})
	if err != nil {
		return err
	}
	right1, err := chart.AddYAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Right,
		ScaleLabel: &chartjs.ScaleLabel{LabelString: "soft-clips", Display: chartjs.True}})
	if err != nil {
		return err
	}
	left2, err := chart.AddYAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Left,
		ScaleLabel: &chartjs.ScaleLabel{LabelString: "insert-size", Display: chartjs.True}})
	if err != nil {
		return err
	}
	left1, err := chart.AddYAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Left,
		ScaleLabel: &chartjs.ScaleLabel{LabelString: "depth", Display: chartjs.True}})
	if err != nil {
		return err
	}

	if _, err = chart.AddXAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Bottom, Display: chartjs.True, ScaleLabel: &chartjs.ScaleLabel{LabelString: "genomic position", Display: chartjs.True}}); err != nil {
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

	gz := gzip.NewWriter(w)
	buf, err := json.Marshal(chart)
	if err != nil {
		return err
	}
	_, err = gz.Write(buf)
	if err != nil {
		return err
	}
	gz.Close()

	return nil
}

func main() {

	cli := &cliarg{}
	cli.Options.MinBaseQuality = 10
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
