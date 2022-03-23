awk 'BEGIN{ chr=$1;\
            pos=$2;\
          }\
          {\
            cur_chr=$1
            cur_pos=$2
            if (cur_chr==chr) {\
              if(cur_pos<pos){\
                ;\
              }\
              else {\
                print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$6;\
              }\
            }\
            pos=cur_pos;\
            chr=cur_chr;\
          }\
          END {\
           print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 ;\
          }' $1
