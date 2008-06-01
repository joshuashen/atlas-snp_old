# input format: 
# read_ID    chromosome    start    end    direction

require 'getoptlong'

opts = GetoptLong.new(
    ["--map", "-m", GetoptLong::REQUIRED_ARGUMENT],
    ["--gsize", "-s", GetoptLong::OPTIONAL_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

if optHash.key?("--help") 
  $stderr.puts "\nUsage: ruby __.rb -m reads.map [ -s genome_size] > output \n\n This script must be executed on single run"
  $stderr.puts "  -m  reads.map, parsed mapping file default format: \n read_ID    chromosome    start    end    direction"
  $stderr.puts "  -s  optional, the size of target genome. It's useful to estimate the expected number of duplicated reads by chance."

  $stderr.puts " \n*Estimating the expected number of reads with same start positions*"
  $stderr.puts " Let N be the number of reads in one run, s be the genome size, 
imagine the start positions of the reads, then based on Lander-Waterman model,
 the expected number of 'island' starts is s*exp(-N/s), therefore the number 
of reads with same start positions is N - (s - s*exp(-N/s)) . A typical 454 
run produces 480000 reads. For human, this number is much smaller than genome
size, so the expected number of same-start reads is close to zero. On the other 
hand, for microbes, such as E. coli, the expected number is about 27000.\n\n"
  exit
end

$hashMaker = proc {|h,k| h[k]=Hash.new(&$hashMaker)} # unlimited-dimensional hash

$hits =Hash.new(&$hashMaker) 

num = 0
File.new(optHash["--map"], "r").each do |line|
  cols = line.split(/\s+/)
  read, chr, s,e,dir = cols[0], cols[1], cols[2].to_i, cols[3].to_i, cols[4]
  
  if dir == '-' or dir == '-1' or dir == 'R' or dir == 'r'
    $hits[chr][e][:r][read] = s
  else
    $hits[chr][s][:f][read] = e
  end
  num += 1
end

fdup = 0
rdup = 0
$hits.each_key do |chr|
  $hits[chr].keys.sort.each do |s|
    if $hits[chr][s][:f].keys.size > 1 
      puts "#{chr}\t+\t#{s}\t#{$hits[chr][s][:f].keys.size}\t#{$hits[chr][s][:f].keys.join(",")}"
      fdup += $hits[chr][s][:f].keys.size - 1
    end
    if $hits[chr][s][:r].keys.size > 1
      puts "#{chr}\t-\t#{s}\t#{$hits[chr][s][:r].keys.size}\t#{$hits[chr][s][:r].keys.join(",")}"
      rdup += $hits[chr][s][:r].keys.size - 1
    end

  end
end

$stderr.puts "Number of duplicated reads in + strand: #{fdup}"
$stderr.puts "Number of duplicated reads in - strand: #{rdup}"
$stderr.puts "Number of duplicated reads in total : #{fdup + rdup}"

if optHash.key?("--gsize") 
  gsize = optHash["--gsize"].to_f
  e = (num - gsize + gsize*Math.exp(-1 * num / gsize)).round
  $stderr.puts "Number of same-start reads expected by chance: #{e}"
else
  
  $stderr.puts "Number of same-start reads expected by chance: #{num} - s + s*exp(-1*#{num}/s)"
  $stderr.puts " s:  the size of the target genome in bp."
end


## In 454 sequencing, some shotgun fragments are occasionally duplicated in the emulation PCR step. This creates skewed data for diploid genomes at those heterozygous loci. To avoid this one needs to remove the duplicated copies from SNP calling. 
## The cause of duplicated reads (from George): In production, When you make the sheared fragments in the first place, you only end up with a limited number of fragments. Then when you do emulsion PCR (for 454), or the amplification to make clusters for Solexa, sometimes fragments will escape from the micelle and populate new beads. But the new bead now has a duplicate fragment. So when you sequence, if you sequence "too much" from a library (imagine sequencing everything) you pick up these duplicates. Our solution has been to make multiple libraries and sequence less from each one. That seems to help a lot.

# So what one would like is to QC libraries for duplicates to make sure you haven't sequenced too much or that there hasn't been too much escape from micelles to form duplicates.


# In practice, duplicated reads are identified by their  start positions on the reference. In a run (the unit of a sequencing experiment) with R reads from a target sequence with size N, the expected number of duplicated reads by pure chance is:  R-N+N*exp(-R/N)






