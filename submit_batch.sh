mkdir -p output
./ditto_xargs jobs.txt
hadd -f output.root output/*root
