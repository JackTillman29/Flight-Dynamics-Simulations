BEGIN{
icol=get_column
}

//{

# separate header into array
nitems = split($0,header_items,",")
printf("%s\n",header_items[icol])

# stop execution
exit

}
