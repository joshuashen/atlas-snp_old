# input: cross_match -discrep_lists output and the reference sequence
# output:  raw SNP list with relevant information

# Author: Yufeng "Joshua" Shen
# Dec, 2008
# 
# What's new: 
#    (perhaps) less memory usage


require 'getoptlong'

# global variables
$basechange={"A"=>"T","T"=>"A","G"=>"C","C"=>"G","N"=>"N","X"=>"X","*"=>"*"}
# cross_match result pattern
$pattern=/^\s*(\d+)\s+(\d+\.\d+)\s+([0-9\.]+)\s+([0-9\.]+)\s+(\S+)\s+(\d+)\s+(\d+)\s+\((\d+)\)\s+(C*?)\s*?(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\**?)/

$coverage = {}
$snp = {}
$seq = {}
$options = {:maxsub => 5.0, :maxindel => 5.0, :minscore => 30} # default options

def main
  optHash = getArguments()
  if optHash.key?("--help") 
    help()
    exit
  end
  setOptions(optHash)
  initiate(optHash)
  scan(optHash)
  output(optHash)  
  if optHash.key?("--coverage")
    coverage(optHash)
  end
end

def getArguments() 
  opts = GetoptLong.new(
                        ["--crossmatch", "-x", GetoptLong::REQUIRED_ARGUMENT],
                        ["--reference", "-r", GetoptLong::REQUIRED_ARGUMENT],
                        ["--output", "-o", GetoptLong::REQUIRED_ARGUMENT],
                        ["--coverage", "-c", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--minscore", "-m", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--maxsub", "-s", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--maxindel", "-g", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  return optHash
end

# $stderr.puts "Note: This program requires large RAM usage. Please run it on servers equipped with large memory."
def help()
  $stderr.puts "Usage: ruby __.rb -x cross_match_output -r reference.fasta -o prefix_of_output [-m minscore -s max_substitution -g max_gap] [-c]"
end

def setOptions(optHash)
  if optHash.key?("--maxsub")
    $options[:maxsub] = optHash["--maxsub"].to_f
  end
  
  if optHash.key?("--maxindel")
    $options[:maxindel] = optHash["--maxindel"].to_f
  end
  
  if optHash.key?("--minscore")
    $options[:minscore] = optHash["--minscore"].to_f
  end
end

def initiate(optHash)
  ref = ''
  File.new(optHash["--reference"], 'r').each do |line|
    if line=~ /^>(\S+)/
      ref = $1
      $seq[ref] = ''
      $snp[ref] = {}
    else
      $seq[ref] << line.chomp
    end
  end
  
  # initiate $coverage
  $seq.each_key do |ref|
    #  $stderr.puts ref
    $coverage[ref] = Array.new($seq[ref].size + 1) {|i|  0}
  end
end

  
def compute(name, ref, span, snps)
  return if span.length < 1

  span.sort! {|a,b| a[0] <=> b[0]}
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
    #    $stderr.puts "#{ref}\t#{s}\t#{e}"
    (s..e).each do |i|
      $coverage[ref][i] += 1
    end
  end

  snps.each_key do |pos|
    refbase = $seq[ref][pos-1,1].upcase
    curbase = snps[pos][:snpbase]
    if snps.key?(pos + 1) or snps.key?(pos + 2) or snps.key?(pos - 1) or snps.key?(pos - 2)
      if snps.key?(pos + 1) and refbase == snps[pos+1][:snpbase] and curbase == $seq[ref][pos,1].upcase
        snps[pos][:info] << "swap;"
      elsif snps.key?(pos - 1) and refbase == snps[pos-1][:snpbase] and curbase == $seq[ref][pos-2,1].upcase 
        snps[pos][:info] << "swap;"
      elsif snps.key?(pos + 2) and refbase == snps[pos+2][:snpbase] and curbase == $seq[ref][pos+1,1].upcase
        snps[pos][:info] << "swap;"
      elsif snps.key?(pos - 2) and refbase == snps[pos-2][:snpbase] and curbase == $seq[ref][pos-3,1].upcase
        snps[pos][:info] << "swap;"
      elsif snps.key?(pos + 1) or snps.key?(pos - 1) 
        snps[pos][:info] << "mnp;"
      else
        snps[pos][:info] << "snp;"        
      end
    else
      snps[pos][:info] << "snp;"      
    end
    $snp[ref][pos] = '' unless $snp[ref].key?(pos)
    $snp[ref][pos] << snps[pos][:info]
  end
end

# count the occurence of homo-polymer runs in a short sequence
def homocount(str)
  s = str.size
  lasts = ''
  homo = {}
  st = 0
  flag = 0
  0.upto(s-1) do |i|
    cur = str[i,1].upcase
    if cur == lasts
      if flag == 0 # this is the second 
        homo[i-1] = 2
        st = i - 1
      else # this is the third or more
        homo[st] += 1
      end
      flag = 1
    else
      flag = 0
    end
    lasts = cur
  end
  if homo.size > 0
    top = homo.keys.sort {|a,b| homo[b] <=> homo[a]}[0]
    xx = homo[top]
  else
    xx = 1
  end
  return xx
end

def scan(optHash)
  offset = 0
  query,ref,qsize,score,dir = '','', 0, 0,''
  flag = 0
  span = []
  snps = {}
  sub, gap, tail = 0, 0, 0 
  File.new(optHash["--crossmatch"],'r').each do |line|
    if line.match($pattern)
      compute(query,ref,span,snps)
      score,sub,del,ins,query,qstart,qend,qright,tstrand,target,d1,d2,d3,lab =$1.to_f,$2.to_f,$3.to_f,$4.to_f,$5,$6.to_i,$7.to_i,$8.to_i,$9,$10,$11,$12,$13,$14
      
      span  = []
      snps = {}
      
      ref = ''
      qsize = qend + qright
      dir = '+'
      ratio = (qend - qstart + 1)/qsize.to_f
      gap = del + ins
      tail = [qstart, qright].max
      if sub > $options[:maxsub] or gap > $options[:maxindel]
        flag = 0
        next
      else
        flag = 1
      end
      if target=~ /^(\S+)\_(\d+)\_(\d+)$/
        ref = $1
        offset = $2.to_i - 1
      elsif target=~ /^(\S+)$/
        ref = $1
        offset = 0
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

      elsif type =~ /^S/ # substitution
        tstart = tplace + offset
        next if ii == 'N'
        dist = qsize - qplace 
        snps[tstart] = {}
        snps[tstart][:snpbase] = ii
        snps[tstart][:info] = "#{ii}(#{qual})#{query}(#{dist})(#{score}/#{qsize})#{dir}#{env}(#{sub}/#{gap}/#{tail})"
      end
    end
  end
  compute(query,ref,span, snps)
end

def output(optHash)
  snpout = File.new(optHash["--output"]+".SNP.list", 'w')
  snpout.puts "refName\tcoordinate\trefBase\thomopolymer\trefEnv\tcoverage\tSNPBase\tadjustedQual\toriQual\tnumSNPReads\tnumAlterReads\treads_info"
  $snp.keys.sort.each do |ref|
    $snp[ref].keys.sort.each do |pos|
      bases = {}
      t = 0
      $snp[ref][pos].split(';').each do |r|
        if r=~ /^(\S)\((\d+)\)(\S+)\((\d+)\)\((\S+)\/(\d+)\)/
          base = $1
          rname = $3
          t += 1
          if !bases.key?(base)
            bases[base] = {}
            bases[base][:adjQual] = 0
            bases[base][:oriQual] = 0
            bases[base][:num] = 0
          end
          bases[base][:num] += 1
          phredQual = $2.to_i
          bases[base][:oriQual] += phredQual
          dist = $4.to_i
          qsize = $6.to_i
          
          
          adjqual1 = phredQual
          bases[base][:adjQual] += adjqual1
        end
      end
      array = bases.keys.sort {|a,b| bases[b][:num] <=> bases[a][:num]}
      #  alternative  = array.size - 1
      qual = bases[array[0]][:adjQual]
      oriqual = bases[array[0]][:oriQual]
      num = bases[array[0]][:num]
      refBase = $seq[ref][pos-1,1]
      refEnv = $seq[ref][pos-7,13]
      homonum = homocount(refEnv)
      
      snpout.puts "#{ref}\t#{pos}\t#{refBase}\t#{homonum}\t#{refEnv}\t#{$coverage[ref][pos]}\t#{array[0]}\t#{qual}\t#{oriqual}\t#{num}\t#{t}\t#{$snp[ref][pos]}"
    end
  end
  snpout.close
end

def coverage(optHash)
  covout = File.new(optHash["--output"]+".ref.depth_cov", 'w')

  $coverage.keys.sort.each do |ref|
    covout.syswrite(">#{ref} depth-coverage\n")
    s = $coverage[ref].size 
    0.upto(s/80) do |i|
      covout.syswrite($coverage[ref][(i*80+1)..(i+1)*80].join(' ') + "\n")
    end
  end
  covout.close
end

main()
