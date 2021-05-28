RiboQC.R \
  --reads ${readFeather} \
  --annotations ${annotationFeather} \
  --whitelist ${whitelist}

rm -f Rplots.pdf
