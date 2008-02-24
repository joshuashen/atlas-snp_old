# compute coverage for targeted region:
#  if the fasta is given, compute the coverage for each base position of the fasta
#  if chromosome and s/e of regions are given, compute the targeted regions


require 'getoptlong'

opts = GetoptLong.new(
    ["--crossmatch", "-x", GetoptLong::REQUIRED_ARGUMENT],
    ["--output", "-o", GetoptLong::OPTIONAL_ARGUMENT],
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
  $stderr.puts "Usage: ruby __.rb -x cross_match.output  -o prefix_of_output [-r ref.fasta] [ -t targeted_genomic_regions ] [ options ]"
  $stderr.puts "    Note: either -r or -t is required "
  $stderr.puts "    Options: "
  $stderr.puts "          -m  -s -g "
  exit
end

if optHash.key?("--output")
  $prefix = optHash["--output"]
else
  $prefix =  optHash["--crossmatch"]   + "_Coverage"
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

$feature = Hash.new {|h,k| h[k] = Hash.new }

def initCov(ref, s, e)
  $coverage[ref] = {} unless $coverage.key?(ref)
 #  $stderr.puts "#{ref}\t#{s}\t#{e}"
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
    if line=~ /^(\S+)\s+(\S+)\s+(\d+)\s+(\d+)/
      featurename, ref,s,e = $1,$2, $3.to_i, $4.to_i
      $feature[ref][featurename] = [s,e]
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



refout = File.new($prefix + '.reference_base_cov', "w")


$coverage.keys.sort.each do |ref|
##  refout.puts ">#{ref}"
  $coverage[ref].keys.sort {|a,b| a<=>b}.each do |i|
    refout.puts "#{ref}\t#{i}\t#{$coverage[ref][i]}"
  end

end

refout.close

refsum = File.new($prefix + ".reference_summary", "w")
$coverage.each_key do |ref|
  tb = $coverage[ref].size
  cc = $coverage[ref].values.select {|i| i>0 }.size

  refsum.puts "#{ref}\t#{tb}\t#{cc}"
end

refsum.close

## statistics

module Math
  module Statistics
    def self.append_features(mod)
      unless mod < Enumerable
        raise TypeError, 
        "`#{self}' can't be included non Enumerable (#{mod})"
      end
      
      def mod.default_block= (blk)
        self.const_set("STAT_BLOCK", blk)
      end
      
      def mod.default_block
        defined?(self::STAT_BLOCK) && self::STAT_BLOCK
      end
      
      super
    end
    
    def default_block
      @stat_block || self.class.default_block
    end
    
    def default_block=(blk)
      @stat_block = blk
    end
    
    def sum
      sum = 0.0
      if block_given?
        each{|i| sum += yield(i)}
      elsif default_block
        each{|i| sum += default_block[*i]}
      else
        each{|i| sum += i}
      end
      sum
    end
    
    def average(&blk)
      sum(&blk)/size
    end
    
    
    def variance(&blk)
      sum2 = if block_given?
               sum{|i| j=yield(i); j*j}
             elsif default_block
               sum{|i| j=default_block[*i]; j*j}
             else
               sum{|i| i**2}
             end
      sum2/size - average(&blk)**2
    end
    
    def standard_deviation(&blk)
      Math::sqrt(variance(&blk))
    end

    
    def Min(&blk)
      if block_given?
        if min = find{|i| i}
          min = yield(min)
          each{|i|
            j = yield(i)
            min = j if min > j
          }
          min
        end
      elsif default_block
        if min = find{|i| i}
          min = default_block[*min]
          each{|i|
            j = default_block[*i]
            min = j if min > j
          }
          min
        end
      else
        min()
      end
    end
    
    def Max(&blk)
      if block_given?
        if max = find{|i| i}
          max = yield(max)
          each{|i|
            j = yield(i)
            max = j if max < j
          }
          max
        end
      elsif default_block
        if max = find{|i| i}
          max = default_block[*max]
          each{|i|
            j = default_block[*i]
            max = j if max > j
          }
          max
        end
      else
        max()
      end
    end
    
    alias avg average
    alias std standard_deviation
    alias var variance
  end
end


class Array
  include Math::Statistics
end


if optHash.key?("--target")
  tout = File.new($prefix + ".target_summary", "w")

  $feature.each_key do |ref|
    $feature[ref].each_key do |ff|
      s, e = $feature[ref][ff]

      cova = $coverage[ref].reject {|p, c| p <  s-1 or p > e - 1}.values   
      avg = cova.average
      stderr = cova.std
      tout.puts "#{ff}\t#{ref}\t#{s}\t#{e}\t#{avg}\t#{stderr}"
    end
  end
  tout.close
end


