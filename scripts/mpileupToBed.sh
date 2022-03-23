threshold=$2;

awk -v minD=$threshold 'BEGIN{ first_start=$2;\
                                final_end=first_start;\
                                coverage=$4;\
                                coveredLength=0;\
                                minCover=coverage;\
                       }\
                       {\
                            current_pos=$2;\
                            coverage=$4;\
                            current_chr=$1;\
                            if(chr!=current_chr){\
                               if(minCover > minD) print chr"\t"first_start"\t"final_end"\t"minCover"\t"coveredLength+1;\
                               first_start=current_pos;\
                               final_end=first_start;\
                               coveredLength=0;\
                               minCover=0;\
                            }\
                            else {\
                              if (current_pos > (final_end + 1)) {\
                               if(minCover > minD) print chr"\t"first_start"\t"final_end"\t"minCover"\t"coveredLength+1;\
                               first_start=current_pos;\
                               final_end=first_start;\
                               coveredLength=0;\
                               minCover=0;\
                             }\
                             else {\
                               final_end=current_pos;\
                               coveredLength=coveredLength+1;\
                               if(minCover < coverage) {\
                                 minCover=coverage;\
                               }\
                             }\
                           }\
                           chr=current_chr;\
                       }\
                       END {\
                       print chr"\t"first_start"\t"final_end"\t"minCover"\t"coveredLength+1;\
                       }' $1
           

          

