##  convert into a format that can be used by R to do glm or support vector machine

# features:
#  1. homopolymer size
#  2. original quality score
#  3. distance to 3' end
#  4. basechange type

# todo:  
#  upstream or downstream of homopolymer? 
##  base swap (consecutive or 1bp apart)? 



trueS = ARGV[0]

falseS = ARGV[1]

puts "class\tswap\thomo\tqual\tdist3\ttail\tsub\tindel\tbasechange"

File.new(trueS, "r").each do |line|
  cols = line.split(/\s+/)

  swap, basechange, homo, qual, dist3, tail, sub, indel, size =  cols[1].to_f, cols[2].to_i, cols[3].to_i, cols[5].to_i, cols[6].to_i,cols[7].to_f, cols[8].to_f,cols[9].to_f, cols[11].to_i

#   dist5 = size - dist3


  oridist = dist3
  if dist3 > 120
    dist3 = 120
  end

  puts "1\t#{swap}\t#{homo}\t#{qual}\t#{dist3}\t#{tail}\t#{sub}\t#{indel}\t#{basechange}"
end


File.new(falseS, "r").each do |line|
  cols = line.split(/\s+/)
  swap, basechange, homo, qual, dist3, tail, sub, indel, size =  cols[1].to_f, cols[2].to_i, cols[3].to_i, cols[5].to_i, cols[6].to_i,cols[7].to_f, cols[8].to_f,cols[9].to_f, cols[11].to_i

  
#  dist5 = size - dist3

  oridist = dist3
  if dist3 > 120
    dist3 = 120
  end
  puts "0\t#{swap}\t#{homo}\t#{qual}\t#{dist3}\t#{tail}\t#{sub}\t#{indel}\t#{basechange}"

end


