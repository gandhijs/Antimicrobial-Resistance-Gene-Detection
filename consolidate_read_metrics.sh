METRICS_DIR=/home/gandhijs/isolates_flagged_mcr-1/CG_pipeline_results

if [ ! -f $METRICS_DIR/master_read_metrics.csv ]
    then touch $METRICS_DIR/master_read_metrics.csv
    echo "File,avgReadLength,totalBases,minReadLength,maxReadLength,avgQuality,numReads,PE?,coverage,readScore,medianFragmentLength" >> $METRICS_DIR/master_read_metrics.csv
fi

for metrics_file in $METRICS_DIR/*.txt
do
    # writes the second line from the read metrics file to a master csv file
    cat $metrics_file | awk 'NR==2' | tr "\\t" "," >> $METRICS_DIR/master_read_metrics.csv 
    
done
