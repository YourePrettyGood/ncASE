#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
}
!/^#/{
#Skip indels and >2 alleles, only consider biallelic SNPs and invariant sites:
   if (length($4) == 1 && length($5) == 1) {
#Identify the element of the FORMAT column corresponding to GT:
      split($9, formatarr, ":");
      for (formatelem in formatarr) {
         if (formatarr[formatelem] == "GT") {
            gtelem = formatelem;
            break;
         };
      };
      split($10, samplearr, ":");
      gt = samplearr[gtelem];
#Select sites without missing genotypes (excluding homref for some reason):
      if (gt == "0/1" || gt == "1/1") {
#Extract their DP4 so that we can calculate total support for each allele:
         split($8, infoarr, ";");
         for (infoelem in infoarr) {
            split(infoarr[infoelem], tagarr, "=");
            if (tagarr[1] == "DP4") {
               split(tagarr[2], dpfour, ",");
               print $1, $2, $4, $5, dpfour[1] + dpfour[2], dpfour[3]+ dpfour[4];
            };
         };
      };
   };
}
