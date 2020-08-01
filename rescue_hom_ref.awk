#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
}
#First file is TSV of union of AIMs:
FNR==NR{
   sites[$1,$2]=toupper($3) SUBSEP toupper($4);
}
#Second file is VCF
FNR<NR&&!/^#/{
#Skip sites not in the AIMs list:
   if (($1, $2) in sites) {
      
   } else {
      next;
   };
#Skip indels and >2 alleles, only consider biallelic SNPs and invariant sites:
#Maybe later we can consider recovering info from multiallelic SNPs
# when the union of AIMs includes such a multiallelic site?
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
#Select sites without missing genotypes:
      if (gt == "0/0" || gt == "0/1" || gt == "1/1") {
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
         split(sites[$1,$2], alleles, SUBSEP);
         #If alleles match between AIMs file and VCF (or are flipped in
         # VCF), output in VCF order:
         if ((toupper($4) == alleles[1] && toupper($5) == alleles[2]) || (toupper($4) == alleles[2] && toupper($5) == alleles[1])) {
            print $1, $2, toupper($4), toupper($5), refcount, altcount;
         #If alt is missing (i.e. GT==0/0), fill in the missing allele:
         } else if ($5 == "." && toupper($4) == alleles[1]) {
            print $1, $2, toupper($4), alleles[2], refcount, altcount;
         } else if ($5 == "." && toupper($4) == alleles[2]) {
            print $1, $2, toupper($4), alleles[1], refcount, altcount;
         #Error out if the alleles don't match the AIMs file:
         } else {
            print "Allele mismatch at site "$1":"$2 > "/dev/stderr";
            print "VCF: "toupper($4)","toupper($5)" AIMs: "alleles[1]","alleles[2] > "/dev/stderr";
            print "Omitting site" > "/dev/stderr";
            next;
         };
      };
   };
}
