# HDAC4 RNA-seq project
## _Frank Zappulla_

[![N|Solid](https://ddfoqzqsu0zvp.cloudfront.net/media/documents/school-of-engineering-wordmark-side-blue.png)](https://www.cse.uconn.edu/)
---
Workflow works in two continguous parts. Run `workflow.sh` to quantify your transcripts using [Salmon](https://combine-lab.github.io/salmon/). Then the newly-generated counts file is the input to the `workflow.py` script which will filter genes so that only genes that have a minimum count of 10 from any 3 of 10 samples will be kept. Then the counts are normalized to the 75th percentile. Finally log2ratio is calculated and written to a file. 


### To Run Workflow Shell Script

- Adjust **config.sh** file to your environment
- Ensure you have all the tools installed that are named in **config.sh**
- Ensure you have sample files in a directory that correspond to the sample names in **fastq_gz_file_pairs.txt**

```
bash workflow.sh
```

### To Run Workflow Python Script

- Adjust **config.py** file to suit your environment
- Ensure you have all the python packages that are imported at the top of the **workflow.py** file
- If you need to recreate **gene_conversion.json**, then run the **gene_conversion.py** file

```
python workflow.py
```


### License

MIT

**Free Software, Hell Yeah!**
