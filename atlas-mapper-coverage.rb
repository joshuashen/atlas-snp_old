# compute coverage for targeted region:
#  if the fasta is given, compute the coverage for each base position of the fasta
#  if chromosome and s/e of regions are given, compute the targeted regions


require 'getoptlong'

opts = GetoptLong.new(
    ["--crossmatch", "-x", GetoptLong::REQUIRED_ARGUMENT],
    ["--reference", "-r", GetoptLong::OPTIONAL_ARGUMENT],
    ["--target", "-t", GetoptLong::OPTIONAL_ARGUMENT],
    ["--minscore", "-m", GetoptLong::OPTIONAL_ARGUMENT],
    ["--maxsub", "-s", GetoptLong::OPTIONAL_ARGUMENT],
    ["--maxindel", "-g", GetoptLong::OPTIONAL_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT]
)


optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

if optHash.key?("--help") or (!optHash.key?("--reference") && !optHash.key?("--target") )
  $stderr.puts "Usage: ruby __.rb -x cross_match.output [-r ref.fasta] [ -t targeted_genomic_regions ] [ options ]"
  $stderr.puts "    Note: either -r or -t is required "
  $stderr.puts "    Options: "
  $stderr.puts "          -m  -s -g "
  exit
end

if optHash.key?("--maxsub")
  $maxsub = optHash["--maxsub"].to_f
else
  $maxsub = 5.0
end

if optHash.key?("--maxindel")
  $maxindel = optHash["--maxindel"].to_f
else
  $maxindel = 5.0
end

if optHash.key?("--minscore")
  $minscore = optHash["--minscore"].to_f
else
  $minscore = 24.0
end

if optHash.key?("--slope")
  $slope = optHash["--slope"].to_f
else
  $slope = 0.13
end


$pattern=/^\s*(\d+)\s+(\d+\.\d+)\s+([0-9\.]+)\s+([0-9\.]+)\s+(\S+)\s+(\d+)\s+(\d+)\s+\((\d+)\)\s+(C*?)\s*?(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\**?)/

$coverage = {}

$seq = {}


def initCov(ref, s, e)
  $coverage[ref] = {} unless $coverage.key?(ref)
  
  (s..e).each do |i|
    $coverage[ref][i] = 0
  end
end

# compute coverage
def compute(name, ref, span)
  return if span.length < 1
  return unless $coverage.key?(ref)

  head = span.shift
  ss,ee = head[0],head[1]
  array = []
  
  array << ss
  span.each do |breaks|
    array << breaks[0]
    array << breaks[1]
  end
  array << ee
  
  while array.size > 0
    s = array.shift
    e = array.shift
    (s..e).each do |i|
      next unless $coverage[ref].key?(i)
      $coverage[ref][i] += 1
    end
  end
end

if optHash.key?("--target")  # just compute the coverage of target regions
  File.new(optHash["--target"], "r").each do |line|
    if line=~ /^(\S+)\s+(\d+)\s+(\d+)/
      ref,s,e = $1,$2.to_i, $3.to_i
      initCov(ref,s,e)
    end
  end
elsif optHash.key?("--reference")
  ref = ''
  s,e = 0, 0 
  seq = ''
  File.new(optHash["--reference"], 'r').each do |line|
    if line=~ /^>(\S+)/
      if seq != ''
        initCov(ref, 1, seq.size)
      end
      ref = $1
      seq = ''
    else
      seq << line.chomp
    end
  end
  if seq != ''
    initCov(ref, 1, seq.size)
  end
end
seq = nil

query,ref,span = '','', []
flag = 0
offset = 0
File.new(optHash["--crossmatch"],'r').each do |line|
  if line.match($pattern)
    compute(query,ref,span)
    score,sub,del,ins,query,qstart,qend,qright,tstrand,target,d1,d2,d3,lab =$1.to_f,$2.to_f,$3.to_f,$4.to_f,$5,$6.to_i,$7.to_i,$8.to_i,$9,$10,$11,$12,$13,$14

    span  = []
    ref = ''
    qsize = qend + qright
    dir = '+'
    ratio = (qend - qstart + 1)/qsize.to_f

    if sub > $maxsub or del + ins > $maxindel
      flag = 0
      next
    else
      flag = 1
    end
    if target=~ /^(\S+)\_(\d+)\_(\d+)$/
      ref = $1
      offset = $2.to_i - 1
    end
    
    if d1 =~ /\((\d+)\)/
      tstart, tend = d2.to_i, d3.to_i
    elsif d3=~ /\((\d+)\)/
      tstart, tend  = d1.to_i, d2.to_i
    end
    if tstart > tend
      tstart, tend = tend, tstart  # wicked                                                                          
      dir = '-'
    end
    s = tstart + offset
    e = tend + offset
    
    #    compute(query,s,e) 
    span << [s,e]
    
  elsif flag > 0 and line=~ /^\s(\S+)\s+(\d+)\s+(\S+)\((\d+)\)\s+(\d+)\s+(\S+)/
    type, qplace, tplace = $1, $2.to_i, $5.to_i
    ii = $3
    qual = $4
    env = $6
    
    if dir == '-'
      ii = $basechange[ii]
    end
    if type =~ /^D/
      # deletion on read                                                                                             
      block = 1
      if type =~ /^D\-(\d+)/
        block = $1.to_i
      end
      
      tstart = tplace
      if dir == '-'  #                                                                             
        tstart = tstart - block + 1
      end
      tend = tstart + block
      span << [tstart + offset, tend + offset]
    end
  end
end
compute(query,ref,span)

$coverage.keys.sort.each do |ref|
  puts ">#{ref}"
 # $coverage[ref].shift  # remove the first one                                                       
  # t = $coverage[ref].size / 5000                                                                                    
 # str = ''                                                                                                          
 # 0.upto(t) do |i|                                                                                                  
    #           covout.puts $coverage[ref][i*5000, 5000].join(" ")                                                   
  # str << "#{$coverage[ref].slice!(i*5000, 5000).join(" ")}\n"                                                    
 # end                                                                                                               
 # covout.puts str          
  $coverage[ref].keys.sort {|a,b| a<=>b}.each do |i|
    puts "#{i}\t#{$coverage[ref][i]}"
    
  end
  #  $coverage[ref]=nil
end

# average depth-cov and completeness

$coverage.each_key do |ref|
  tb = $coverage[ref].size
  cc = $coverage[ref].values.select {|i| i>0 }.size
  
  $stderr.puts "#{ref}\t#{tb}\t#{cc}"
end

