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
         adf="";
         adr="";
         dpfour="";
         split($8, infoarr, ";");
         for (infoelem in infoarr) {
            split(infoarr[infoelem], tagarr, "=");
            if (tagarr[1] == "DP4") {
               dpfour=tagarr[2];
            } else if (tagarr[1] == "ADF") {
               adf=tagarr[2];
            } else if (tagarr[1] == "ADR") {
               adr=tagarr[2];
            };
         };
         #Use DP4 when available:
         if (length(dpfour) > 0) {
            split(dpfour, asd, ",");
            refcount=asd[1] + asd[2];
            altcount=asd[3] + asd[4];
         #Otherwise use ADF and ADR:
         } else if (length(adf) > 0 && length(adr) > 0) {
            refcount=adf;
            altcount=adr;
         #Error out if neither is available:
         } else {
            print "Unable to detect allele depth in VCF at site "$1":"$2", omitting site." > "/dev/stderr";
            next;
         };
         print $1, $2, toupper($4), toupper($5), refcount, altcount;
      };
   };
}
