
while line = ARGF.gets do 
  if line=~ /^\s*S\s+(\d+)\s+\S+\((\d+)\)\s+\d+/
    puts "#{$1}\t#{$2}"
  end
end

