function output_str = strsplit_select1(input_str,delimiter,splitindx)

output_str = strsplit(input_str,delimiter);
output_str = output_str(splitindx);