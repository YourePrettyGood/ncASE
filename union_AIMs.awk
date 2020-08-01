#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
   #Output in a sorted order (not optimal order though):
   PROCINFO["sorted_in"]="@ind_str_asc";
}
#All files are .snps files, so columns 1 and 2 are scaffold and position,
# columns 3 and 4 are ref and alt alleles, and columns 5 and 6 are
# ref and alt allele depths.
{
   if (($1,$2) in aims) {
      split(aims[$1,$2], alleles, SUBSEP);
      if ($3 == alleles[1] && $4 == alleles[2]) { #Same orientation
         
      } else if ($4 == alleles[1] && $3 == alleles[2]) { #Opposite orientation
         
      } else { #Invalid orientation, mismatched alleles
         if (length(debug) > 0) {
            print $1":"$2 " mismatched alleles "alleles[1]","alleles[2]" vs. "$3","$4 > "/dev/stderr";
         };
         #Omit if mismatched?
         if (length(omit_mm) > 0) {
            delete aims[$1,$2];
         };
      };
   } else {
      aims[$1,$2]=$3 SUBSEP $4;
   };
}
END{
   for (i in aims) {
      split(i, pos, SUBSEP);
      split(aims[i], alleles, SUBSEP);
      #Output is *not* sorted normally, could pipe to
      # sort -k1,1V -k2,2n
      #However, that shouldn't be necessary, since
      # rescuing will just use a hash
      print pos[1], pos[2], alleles[1], alleles[2];
   };
}
