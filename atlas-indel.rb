## input: cross_match -discrep_lists output

# limitations: 
## 1.  cross_match one-block alignment only
#  2.  ignoring 1bp indels by default

# the quality of the indel: 
#   the accuracy of the gap size is important, therefore the quality of indels should reflect the confidence of the gap size


require 'getoptlong'

opts = GetoptLong.new(
    ["--crossmatch", "-x", GetoptLong::REQUIRED_ARGUMENT],
    ["--reference", "-r", GetoptLong::REQUIRED_ARGUMENT],
    ["--minscore", "-m", GetoptLong::OPTIONAL_ARGUMENT],
    ["--maxsub", "-s", GetoptLong::OPTIONAL_ARGUMENT],
    ["--maxindel", "-g", GetoptLong::OPTIONAL_ARGUMENT],
    ["--cutoff", "-c", GetoptLong::OPTIONAL_ARGUMENT],
    ["--minIndelSize", "-l", GetoptLong::OPTIONAL_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

if optHash.key?("--help")
  $stderr.puts "Usage: ruby __.rb -x cross_match_result -r reference.fasta [options] > output"
  $stderr.puts "  Warning: it's slow."
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
  $maxindel = 25.0
end

if optHash.key?("--minscore")
  $minscore = optHash["--minscore"].to_f
else
  $minscore = 24.0
end

if optHash.key?("--minIndelSize")
  $minIndelSize = optHash["--minIndelSize"].to_i
else
  $minIndelSize = 2
end

# $now = "atlas-indel-piles-" + Time.now.to_f.to_s + ".txt" 
# $tempfile =File.new($now,"w")

$seq = {}

# for unlimited-dimensional hash
hashMaker = proc {|h,k| h[k]=Hash.new(&hashMaker)}

$deletions = {} # Hash.new(&hashMaker)
$insertions = {} # Hash.new(&hashMaker)

$cov = {}
# $atlas = {} ## stores the mapping information of all reads

class Indel 
  attr_accessor :ref, :query, :tstart, :tend, :qstart, :qend, :dir,:sub, :dels,:ins,:envs, :type, :bases, :dist3
  
  def initialize(query, ref, qstart, qend, tstart, tend, dir,sub,dels,ins,l,envs,type)
    @query,@ref,@qstart,@qend, @tstart, @tend, @dir,@sub,@dels,@ins,@length,@envs,@type = query, ref, qstart, qend, tstart, tend, dir,sub,dels,ins,l,envs,type
    
    @dist3 = @length - @qend
    ts = @tend - @tstart
    qs = @qend - @qstart
    gap = ts

    @bases = ''
    if @type == 'D'
      gap = ts
      @bases = $seq[@ref][@tstart-1..@tend-2]
    else
      gap  = qs
      @bases = @envs[6..-7]
      if @dir == '-'
        @bases = bases.reverse.tr('ATCG', 'TAGC')
      end
    end
  end
  
  def dump
#    puts "#{type}\t#{gap}\t#{@ref}\t#{@tstart}\t#{@tend}\t#{@dir}\t#{@query}\t#{@qstart}\t#{@qend}\t#{@sub}\t#{@dels}\t#{@ins}\t#{@length}\t#{envs}\t#{bases}"
    return "#{@bases},#{@dist3},#{@query},#{@dir},#{@qstart},#{@qend},#{@sub},#{@dels},#{@ins},#{@length},#{@envs}"
  end
end

def assemble(indelarray)
  # goal: get consensus,  

# rules:

  bases = {}

#!! need to sort by distance to 3' end -- the larger the more accurate
#  indelarray.sort {|a,b| a.dist3 <=> b.dist3 }.each do |indel|
  indelarray.each do |indel|
    if bases.key?(indel.bases)
      bases[indel.bases] +=1 
    else
      bases[indel.bases] = 1
    end
  end

  consensus = bases.keys.sort {|a,b| bases[a]<=>bases[b]}[-1]

#  gap = consensus.length

  return consensus
end

# store the coverage information in a temporary file to save RAM
def pileup(span, ref)

  return if span.size < 1 or !$seq.key?(ref)
  # $stderr.puts "#{span.join(' ')}\t#{ref}"
  head = span.shift
  ss, ee = head[0], head[1]
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
  #  $stderr.puts "#{ref}\t#{s}\t#{e}"   
#    $atlas[ref] << [s,e]
#    $tempfile.puts "#{ref}\t#{s}\t#{e}"
    ((s+1)..(e-1)).each do |i|
      $cov[ref][i] += 1
    end
  end
end

ref = ''
# read ref fasta, initialize

File.new(optHash["--reference"], "r").each do |line|
  if line=~ /^>(\S+)/
    $seq[$1] = ''
    ref = $1
    $insertions[ref] = Hash.new(&hashMaker) 
    $deletions[ref] = Hash.new(&hashMaker) 
#    $atlas[ref] = []
  else
    $seq[ref] << line.chomp
  end
end

$seq.each_key do |ref|
#  $stderr.puts ref                                                                                
  $cov[ref] = Array.new($seq[ref].size + 1) {|i|  0}
end


pattern=/^\s*(\d+)\s+(\d+\.\d+)\s+([0-9\.]+)\s+([0-9\.]+)\s+(\S+)\s+(\d+)\s+(\d+)\s+\((\d+)\)\s+(C*?)\s*?(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\**?)/

sub,query, target, dir = 0.0,'','',0
dels,ins = 0.0, 0.0
length = 0
flag = 0
offset = 0
ref = ''
span = []

File.new(optHash["--crossmatch"], "r").each do |line|
  if line.match(pattern)
    pileup(span, ref)
    score,sub,dels,ins,query,qstart,qend,qright,tstrand,target,d1,d2,d3,lab =$1.to_f,$2.to_f,$3.to_f,$4.to_f,$5,$6.to_i,$7.to_i,$8.to_i,$9,$10,$11,$12,$13,$14
    span  = []
    if sub > $maxsub or dels + ins > $maxindel 
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
      tright, tstart, tend = $1.to_i, d2.to_i, d3.to_i
    elsif d3=~ /\((\d+)\)/
      tstart, tend, tright = d1.to_i, d2.to_i, $1.to_i
    end
      
    length = qright + qend
    if tstart < tend
      dir = '+'
    else
      dir = '-'
      tstart, tend = tend, tstart  # wicked 
    end
    span << [tstart, tend]
  elsif flag > 0 and line=~ /^\s(\S+)\s+(\d+)\s+(\S+)\(\d+\)\s+(\d+)\s+(\S+)/
    type, qplace, tplace,envs = $1, $2.to_i, $4.to_i,$5
#    if type =~ /^D\-(\d+)/
    if type =~ /^D/
      # deletion on read
      qstart = qplace 
      qend = qstart 
      block = 1
      if type =~ /^D\-(\d+)/
        block = $1.to_i
      end
      
      next if block < $minIndelSize
      
      tstart = tplace 

      if dir == '-'  #
        tstart = tstart - block + 1
      end
      tend = tstart + block

      tstart = tstart + offset
      tend = tend  + offset
      indel = Indel.new(query,ref, qstart,qend, tstart, tend, dir,sub,dels,ins,length,envs, "D")

      
      $deletions[ref][tstart][indel] = 1
#      $deletions[ref][tstart][:through] = 0
#      $deletion[ref][tstart][:b2] = 0
      span << [tstart, tend]
    elsif type =~ /^I/
      qstart = qplace 
      tstart = tplace  
      block = 1
      if type =~ /^I\-(\d+)/
        block = $1.to_i
      end

      next if block < $minIndelSize

      qend = qstart + block
     
      if dir == '-'
        tstart -= 1
      end
      tend = tstart

      tstart = tstart +offset
      tend = tend  + offset
      
      indel = Indel.new(query,ref, qstart,qend, tstart, tend, dir,sub,dels,ins,length,envs,"I")
      $insertions[ref][tstart][indel] = 1
#      $insertions[ref][tstart][:through]  = 0
      span << [tstart, tend + 1]
    end
  end
end
pileup(span, ref)


#File.new($now, "r").each do |line|
#  cols = line.split(/\s+/)

#  ref,s,e = cols[0], cols[1].to_i, cols[2].to_i
# $atlas.each_key do |ref|
#  $atlas[ref].sort! {|a,b| a[0] <=> b[0]}
  # $atlas[ref].shift

#  ((s+1)..(e-1)).each do |i|
#    if $insertions[ref].key?(i)
#   number of reads that walk through the break points
#      $insertions[ref][i][:through] += 1
#    end
    
#    if $deletions[ref].key?(i) and $deletions[ref][i][:consensus].size + i -1 < e  # go through
#      $deletions[ref][i][:through] += 1
#    end
#  end
# end

$seq = nil

# header
puts "indel_type\tindel_size\tchr\tstart\tend\tindelbases\tnumIndelReads\tnumRefReads\tinfo_indelReads"

$insertions.keys.sort.each do |ref|
  $insertions[ref].keys.sort.each do |s|
    consensus = assemble($insertions[ref][s].keys)

    cov = $cov[ref][s]

    print "I\t#{consensus.size}\t#{ref}\t#{s}\t#{s+1}\t-/#{consensus}\t#{$insertions[ref][s].size}\t#{cov}\t"
    $insertions[ref][s].each_key do |indel|
      print "#{indel.dump};"
    end
    print "\n"
  end
end

$deletions.keys.sort.each do |ref|
  $deletions[ref].keys.sort.each do |s|
    c = assemble($deletions[ref][s].keys)

    e = s + c.size - 1
    m = (s+e)/2
    cov = [$cov[ref][s],$cov[ref][e], $cov[ref][m] ].min

    print "D\t#{c.size}\t#{ref}\t#{s}\t#{e}\t#{c}/-\t#{$deletions[ref][s].size}\t#{cov}\t"
    $deletions[ref][s].each_key do |indel|
      print "#{indel.dump};"
    end
    print "\n"
  end
end

# $tempfile.close

# File.delete($now)
