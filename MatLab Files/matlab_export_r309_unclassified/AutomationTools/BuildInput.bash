#!/bin/bash

# is LibreOffice on this machine?
#     need for "unoconv" which is LO's file format converter

#===================================================
# definitions
#===================================================
tooldir=/home/keeneal/ESAMSDevelopmentBranch/utils/AutomationTools
awkdir=$tooldir/scripts_awk
template_dir=$tooldir
#template_file=Template.inp
new_infile_dir=$tooldir/../../input
mkdir $new_infile_dir

input_filename=GNOMEAutomationGenerationAKeene

# spreadsheet files
#   definitely prefer CSV files
ss_basename=GNOMEAutomationGenerationAKeene
ss_main=${ss_basename}

#===================================================
# convert EXCEL spreadsheet to a CSV file
#   "unoconv" is apart of LibreOffice (Open Office)
#===================================================
ss_ext=xlsx
unoconv -d spreadsheet -f csv ${ss_main}.${ss_ext}

ss_ext=csv
ss_main_csv=${ss_main}.${ss_ext}


#===================================================
# get number of lines in CSV file (excluding header)
#===================================================
nrows=`cat $ss_main_csv | wc -l`; nruns=`echo "$nrows - 1" | bc`
ncols=`cat $ss_main_csv | awk -f $awkdir/get_header_csv.awk | wc -l`
echo "# rows: " $nrows
echo "# cols: " $ncols
echo "# runs: " $nruns

#===================================================
# SETUP INPUT FILES
#===================================================
header=`head -n 1 $ss_main_csv`
for ifile in `seq 2 $nrows`; do
  # get row data
  #    head -n K : prints the first K rows of the file
  #    tail -1   : of the first K lines, print the last one (desired row)
  rowdat=`head -n $ifile $ss_main_csv | tail -1`
  echo $rowdat

  # make the input file if this row's first column ("<Build>") is "1"
  if [ "`echo $rowdat | awk -f $awkdir/get_col_from_row_csv.awk --assign get_column=1`" == "1" ]; then

    # create file using <Data Set Name> and <Unique ID> (cols 2 and 3)
    data_set_name=`echo $rowdat | awk -f $awkdir/get_col_from_row_csv.awk --assign get_column=2`
    unique_id=`echo $rowdat | awk -f $awkdir/get_col_from_row_csv.awk --assign get_column=3`
    new_filename=${data_set_name}_${unique_id}_original.inp
    template_file=`echo $rowdat | awk -f $awkdir/get_col_from_row_csv.awk --assign get_column=4`
    echo $template_file
    cp $template_dir/$template_file  $new_infile_dir/$new_filename
    echo "making  ---  " $new_infile_dir "/" $new_filename

    # go through each header item
    # search/replace "<header item" --> "value"
    for icol in `seq 1 $ncols`; do
      header_str=`echo $header | awk -f $awkdir/get_col_from_row_csv.awk --assign get_column=$icol`
      header_val=`echo $rowdat | awk -f $awkdir/get_col_from_row_csv.awk --assign get_column=$icol`
      sed --in-place -e "s,$header_str,$header_val,g" $new_infile_dir/$new_filename
      #echo $header_str "  " $header_val
    done

  fi
done

echo "   FINISHED!"

#===================================================
# get headers
#===================================================
#for header in `cat $ss_main_csv | awk -f $awkdir/get_header_csv.awk`; do
#for header in `cat $ss_main_csv | awk '//{print $0; exit}' | cut -f@ -d","`; do
#  # get all data associated with this header
#  echo $header
#
#
#done



