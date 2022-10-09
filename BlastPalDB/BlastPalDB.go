package blastpaldb

import (
	"fmt"
	config "gophymmr/Config"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"github.com/linxGnu/grocksdb"
	"strconv"
	"math"
	"encoding/json"
)

type BlastOut struct {
	Target 		int64		`json:"target"`
	Score 		float64		`json:"score"`
	Evalue		float64		`json:"evalue"`
	LogEvalue	float64		`json:"log_value"`
	BlastStart	int64		`json:"start"`
	BlastEnd	int64		`json:"end"`

}

type Triplet[T, U, N any] struct {
	A 	T
	B 	U
	C	N
}

type TripOut = Triplet[string, string, []BlastOut]

var db *grocksdb.DB
var dbpath string

func init() {
	bbto := grocksdb.NewDefaultBlockBasedTableOptions()
	bbto.SetBlockCache(grocksdb.NewLRUCache(3 << 30))

	opts := grocksdb.NewDefaultOptions()
	opts.SetBlockBasedTableFactory(bbto)
	opts.SetCreateIfMissing(true)

	var err error

	db, err = grocksdb.OpenDb(opts, config.DEFAULT_DB_PATH)

	if err != nil {
		log.Fatalln(err)
	} else {
		log.Println("Database loaded")
	}
}

func processGeneName(gene, blastPathBase string) (string, string) {
	blastFileName := fmt.Sprintf("%s.blast", gene)
	blastPath := filepath.Join(blastPathBase, blastFileName)
	result := fmt.Sprintf("%s.done", blastPath)

	return result, blastPath
}

func execCommand(gene_hits map[string]string, dbPath, gene, tmp_path, result, blastPath, prog string, evalueThresh float64, scoreThresh, numThread int) string {
	if stat, err := os.Stat(result); os.IsNotExist(err) || stat.Size() == 0 {
		if _, err := os.Stat(blastPath); os.IsNotExist(err) {
			gene_name := fmt.Sprintf("%s.fa", gene)

			tmpPath := filepath.Join(tmp_path, gene_name)

			f, err := os.Create(tmp_path)

			if err != nil {
				log.Fatal(err)
			}

			for target, sequence := range gene_hits {
				f.Write([]byte(fmt.Sprintf(">%s\n%s", target, sequence)))
			}

			if _, err := os.Stat(blastPath); !os.IsNotExist(err) {
				fBlast, errBlast := os.Open(blastPath)

				if errBlast != nil {
					log.Fatalln(errBlast)
				}

				insert := make([]byte, 0)

				fBlast.Write(insert)
			}

			cmd := exec.Command(prog, fmt.Sprintf("-outfmt '7 qseqid sseqid evalue bitscore qstart qend' -evalue '%f' -threshold '%d' -num_threads '%d' -db '%s' -query '%s' -out '%s'", evalueThresh, scoreThresh, numThread, dbPath, tmpPath, blastPath))

			errCmd := cmd.Run().Error()

			log.Printf("Command ran. Possible errors: %s \n", errCmd)

			os.Remove(tmpPath)
		}

		os.Rename(blastPath, result)
	}

	return result

}


func geneOut(gene, result string, minimumScore, minimumEvalue float64) (map[string][]BlastOut, []TripOut) {
	geneOut := make(map[string][]BlastOut, 0)
	retOut := make([]TripOut, 0)

	thisContent, errContent := os.Open(result)

	if errContent != nil {
		log.Fatal(errContent)
	}

	var contentBytes []byte

	_, errRead := thisContent.Read(contentBytes)

	if errRead != nil {
		log.Fatal(errRead)
	}

	contentStr := string(contentBytes)

	if len(contentStr) != 0 {
		for _, line := range strings.Split(contentStr, "\n") {
			fields := strings.Split(line, "\t")

			if len(fields) == 6 {
				queryId := fields[0]
				subjectId := fields[1]
				evalue := fields[2]
				bitScore := fields[3]
				qStart := fields[4]
				qEnd := fields[5]

				bitScoreFloat, err := strconv.ParseFloat(bitScore, 64)

				if err != nil {
					log.Fatal(err)
				}

				evalueFloat, err := strconv.ParseFloat(evalue, 64)

				if err != nil {
					log.Fatal(err)
				}

				var logValue float64

				if bitScoreFloat >= minimumScore && evalueFloat <= minimumEvalue {
					logValue = math.Log(evalueFloat)
				}

				queryIdSplit := strings.Split(queryId, "_hmmid")
				hmmSearchId := queryIdSplit[1]

				subjectIdInt, err := strconv.ParseInt(subjectId, 10, 64)

				if err != nil {
					log.Fatal(err)
				}

				startInt, err := strconv.ParseInt(qStart, 10, 64)

				if err != nil {
					log.Fatal(err)
				}


				endInt, err := strconv.ParseInt(qEnd, 10, 64)

				if err != nil {
					log.Fatal(err)
				}


				blastOut := BlastOut {		
					Target: subjectIdInt,						
					Score: bitScoreFloat,
					Evalue: evalueFloat,
					LogEvalue: logValue,
					BlastStart: startInt,
					BlastEnd: endInt,
				}

				if _, ok := geneOut[hmmSearchId]; !ok {
					geneOut[hmmSearchId] = make([]BlastOut, 0)
				}

				geneOut[hmmSearchId] = append(geneOut[hmmSearchId], blastOut)
							
			}
		}

		for hmmSearchId, blastOutArray := range geneOut {
			key := fmt.Sprintf("blastfor:%s", hmmSearchId)

			data, err := json.Marshal(blastOutArray)

			if err != nil {
				log.Fatal(err)
			}

			dataStr := string(data)

			retOut = append(retOut, TripOut{key, dataStr, blastOutArray})
		}
	}

	fmt.Printf("Blasted: %s \n", gene)


	return geneOut, retOut

}
