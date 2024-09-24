BEGIN{

}


//{

# separate header into array
nitems = split($0,header_items,",")
for(j = 1; j <= nitems; j++){
  # remove white space " " between the comma separator and the "<" on the header item name
  # leading the string (<var name>)
  gsub(/^ */,"",header_items[j]) # before "<"
  gsub(/ *$/,"",header_items[j]) # after  ">"
  printf("%s\n",header_items[j])
}

# stop execution
exit
}
