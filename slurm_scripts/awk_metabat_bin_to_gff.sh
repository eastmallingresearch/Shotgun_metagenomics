awk '{
  if(index($0,">")){
  header=gensub(/.*\.fa\./,"","g",$0);
  bin=gensub(/>/,"ID=","1",$0);
  if(tot){print header,"METABAT","BINS",1,tot,".","+",".",bin;tot=0}}else{tot=tot+length($0)}
} END {print header,"METABAT","BINS",1,tot,".","+",".",bin}' OFS="\t"
