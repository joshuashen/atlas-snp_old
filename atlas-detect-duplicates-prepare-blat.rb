require 'zlib'

blatResult = ARGV[0]

$reads = {}

def processLine(line)
  cols = line.split(/\s+/)
  return if cols.size < 18
  hit = {}
  
  hit[:match], qgap,tgap,dir, query,qsize,hit[:qstart],hit[:qend],hit[:target],hit[:tsize],hit[:tstart],hit[:tend],hit[:blocks] =cols[0].to_i,cols[5].to_i,cols[7].to_i,cols[8],cols[9], cols[10].to_i,cols[11].to_i+1,cols[12].to_i,cols[13],cols[14].to_i,cols[15].to_i+1,cols[16].to_i,cols[17].to_i
  
  if dir == '+'
    ori = 1
  else
    ori = -1
  end
  
  if ori > 0 
    hit[:tstart] -= hit[:qstart] - 1
  else
    hit[:tend] += hit[:qstart] - 1
  end
  

  if !$reads.key?(query)
    $reads[query] = {}
    $reads[query][:score] =  -100
    $reads[query][:hits] = {}
    $reads[query][:size] = qsize
    $reads[query][:dir] = ori
    $reads[query][:tophit] = nil
  end
  
  s = hit[:match] - (qgap + tgap)/2.0                                                                           
  hit[:score] = s
 

  if s >= $reads[query][:score]

    if s > $reads[query][:score]
      $reads[query][:score] = s
      $reads[query][:tophit] = hit
      $reads[query][:dir] = ori
      $reads[query][:hits] = {}
    end

    $reads[query][:hits][hit] = 1  
  end
end


if blatResult =~ /.gz$/
  fo = Zlib::GzipReader.open(blatResult)
else
  fo = File.new(blatResult, 'r')
end

fo.each do |line|
  if line=~/^\d+/
    processLine(line)
    fo.each do |line|
      processLine(line)
    end
  else
    next
  end
end

fo.close

$stderr.puts "Number of reads processed: #{$reads.size}"
$reads.keys.sort.each do |query|
  score = $reads[query][:score]
  if score <= -100
    $stderr.puts "#{query}\tno good hits\t#{score}"
    next
  end

  hit = $reads[query][:tophit]
  num = $reads[query][:hits].keys.size
  print "#{query}\t#{hit[:target]}\t#{hit[:tstart]}\t#{hit[:tend]}\t#{$reads[query][:dir]}"
  if num > 1
    puts "\trepeat"
  else
    puts "\tunique"
  end
end
