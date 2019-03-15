# Format labels to be more aesthetic
# Slightly modified from this google group discussion (written by Brian Diggs)
# https://groups.google.com/forum/#!topic/ggplot2/a_xhMoQyxZ4

 fancy_scientific <- function(l) {
   # turn in to character string in scientific notation
   l <- format(l, scientific = TRUE)
   # quote the part before the exponent to keep all the digits
   l <- gsub("^(.*)e", "'\\1'e", l)
   # turn the 'e+' into plotmath format
   l <- substr(gsub("e", "10^", l),nchar(l)-3, nchar(l)+2)
   l<-gsub("\\+", "", l)
   # return this as an expression
   parse(text=l)
 }

fancy_scientific2 <- function(l) {
  # turn in to character string in scientific notation
  l<-gsub("\\<0>\\", "zero", l, fixed=TRUE)
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  l<-gsub("\\+", "", l)
  # return this as an expression
  parse(text=l)
}
