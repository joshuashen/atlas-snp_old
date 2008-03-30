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

puts "class\tswap\thomo\trhomo\tqual\tdist3\tdist5\ttail\tsub\tindel\tbasechange"

File.new(trueS, "r").each do |line|
  cols = line.split(/\s+/)

  swap, basechange, homo, rhomo, qual, dist3, tail, sub, indel, size =  cols[1].to_f, cols[2].to_i, cols[3].to_i, cols[4].to_i, cols[6].to_i, cols[7].to_i,cols[8].to_f, cols[9].to_f,cols[10].to_f, cols[12].to_i

  dist5 = size - dist3


  oridist = dist3
  if dist3 > 110
    dist3 = 110
  end
  
#  if dist5< 100
#    dist5 = 100
#  end

  puts "1\t#{swap}\t#{homo}\t#{rhomo}\t#{qual}\t#{dist3}\t#{dist5}\t#{tail}\t#{sub}\t#{indel}\t#{basechange}"
end


File.new(falseS, "r").each do |line|
  cols = line.split(/\s+/)
  swap, basechange, homo, rhomo, qual, dist3, tail, sub, indel, size =  cols[1].to_f, cols[2].to_i, cols[3].to_i, cols[4].to_i, cols[6].to_i, cols[7].to_i,cols[8].to_f, cols[9].to_f,cols[10].to_f, cols[12].to_i

  
  dist5 = size - dist3

  oridist = dist3
  if dist3 > 110
    dist3 = 110
  end
  
 # if dist5 < 100
 #   dist5 = 100
 # end


  puts "0\t#{swap}\t#{homo}\t#{rhomo}\t#{qual}\t#{dist3}\t#{dist5}\t#{tail}\t#{sub}\t#{indel}\t#{basechange}"

end


