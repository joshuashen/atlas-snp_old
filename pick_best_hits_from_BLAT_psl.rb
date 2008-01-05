require 'zlib'

blatResult = ARGV[0]
cut = ARGV[1]

if cut == nil
  $cutoffRatio = 0.995
else
  $cutoffRatio = cut.to_f
end

reads = {}


if blatResult =~ /.gz$/
  fo = Zlib::GzipReader.open(blatResult)
else
  fo = File.new(blatResult, 'r')
end

# File.new(blatResult, 'r').each do |line|

fo.each do |line|
  next unless line=~ /^\d+/
  cols = line.split(/\s+/)

  hit = {}
  
  hit[:match], qgap,tgap,dir, query,qsize,hit[:qstart],hit[:qend],hit[:target],hit[:tsize],hit[:tstart],hit[:tend],hit[:blocks] =cols[0].to_i,cols[5].to_i,cols[7].to_i,cols[8],cols[9], cols[10].to_i,cols[11].to_i+1,cols[12].to_i,cols[13],cols[14].to_i,cols[15].to_i+1,cols[16].to_i,cols[17].to_i
  
  if dir == '+'
    ori = 1
  else
    ori = -1
  end
  if !reads.key?(query)
    reads[query] = {}
    reads[query][:score] =  -100
    reads[query][:hits] = {}
    reads[query][:size] = qsize
    reads[query][:dir] = ori
    reads[query][:tophit] = nil
  end
  
  s = hit[:match] - (qgap + tgap)/2.0                                                                           
  hit[:score] = s
  

  if s >=  reads[query][:score]*$cutoffRatio
    reads[query][:hits][hit] = 1  
    if s > reads[query][:score]
      reads[query][:score] = s
      reads[query][:tophit] = hit
      reads[query][:dir] = ori
    end
  end
end

fo.close

reads.keys.sort.each do |query|
  score = reads[query][:score]
  if score <= -100
    $stderr.puts "#{query}\tno good hits\t#{score}"
    next
  end
  num = 0
  hit = reads[query][:tophit]
#  $stderr.puts "#{query}\t#{score}"
  reads[query][:hits].each_key do |h|
    
    if h[:score] > score*$cutoffRatio
      num += 1
    end
  end
  print "#{query}\t#{reads[query][:size]}\t#{hit[:target]}\t#{hit[:tsize]}\t#{hit[:match]}\t#{hit[:qstart]}\t#{hit[:qend]}\t#{hit[:tstart]}\t#{hit[:tend]}\t#{score}\t#{reads[query][:dir]}\t#{hit[:blocks]}"
  if num > 1
    puts "\trepeat"
  else
    puts "\tunique"
  end
end
