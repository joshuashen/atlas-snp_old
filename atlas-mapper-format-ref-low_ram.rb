# format reference sequence:
# 1. generate ooc file for BLAT
# 2. divide reference sequences into 100K overlapping pieces

## TODO:
# 1. consider circular genomes (microbes)
# 2. [perhaps] split huge genomes, such human genome,into smaller pieces for blat, in order to run BLAT with <=2G RAM.

require 'getoptlong'

opts = GetoptLong.new(
    ["--reference", "-r", GetoptLong::REQUIRED_ARGUMENT],

    ["--length", "-l", GetoptLong::OPTIONAL_ARGUMENT],
    ["--freq","-f", GetoptLong::OPTIONAL_ARGUMENT],
    ["--blat", "-b", GetoptLong::OPTIONAL_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT]
)

package = File.dirname(File.expand_path($0))

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

if optHash.key?("--help")
  $stderr.puts "Usage: ruby __.rb -r reference.fasta [-l length_of_pieces -f cutoff_freq_of_11mer] [-b path_to_blat]"
  $stderr.puts " \n  Note: please make sure the path to blat is set in your shell env. If not, you can provide the path by -b option."
  $stderr.puts "This program builds a directory for BLAT/cross_match between a batch of short sequences, such as 454 reads, and the reference sequence."
  $stderr.puts "    -f [1024] specifies the cutoff of 11-mer frequency deemed as over-represented in the reference so BLAT will ignore it during seeding step."
  exit
end

if optHash.key?("--blat")
  blat = optHash["--blat"]
else
  blat = "blat"
end

oocCut = 1024
kmer = 11
pieceLength = 100000 # default 100Kbp
overlap = 5000


if optHash.key?("--freq")
  oocCut = optHash["--freq"].to_i
end

if optHash.key?("--length")
  pieceLength = optHash["--length"].to_i
end

envDir = optHash["--reference"] + '.Env4mapping'
system("mkdir #{envDir}")
ref = optHash["--reference"]
oocOut = envDir + "/11mer.ooc"

# generate ooc file
#if !File.exist?(File.expand_path(blat))
#	$stderr.puts "Cannot find blat"
#	exit
#end

cmd = "#{blat} #{ref} /dev/null /dev/null -tileSize=#{kmer} -makeOoc=#{oocOut} -repMatch=#{oocCut}"
system(cmd)

# break the reference into 900000000 non-overlapping pieces 

# divisions = envDir + "/divisions.index"

divDir = envDir + "/ref_divisions/"

system("mkdir #{divDir}")

cmd = "perl #{package}/divide-fasta-w-qual.pl -f #{ref} -l 900000000 -p #{divDir}div"
system(cmd)


# break the reference into 100k overlapping pieces and store the information in an index file

indexFile = envDir + "/fragments.index"

$index = {:length=>0, :done=>0, :cut=>pieceLength, :overlap=>overlap, :current => nil}

$index[:basedir] = envDir + "/ref_pieces/"

system("mkdir #{$index[:basedir]}")

$index[:num] = 0
$index[:layer] = 0
$index[:index] = {}
$index[:output] = File.new(indexFile, 'w')
$index[:refpiece] = "" 
def formatDB(name, seq, index)
  return if name == '' or seq.length < 10
  
  if index[:length] + seq.length <= index[:cut]   # add the whole seq
    if index[:current] == nil # no file created
      index[:num] += 1
      index[:length] = 0
      if index[:num] % 1000 == 1  # 1000 per dir
        index[:layer] +=1 
        laydir = index[:basedir] + "#{index[:layer]}" + '/'
        system("mkdir #{laydir}")
      end
      index[:refpiece] = index[:basedir] + "#{index[:layer]}" + '/' + 'piece_' + "#{index[:num]}" + '.fa'
      index[:current] = File.new(index[:refpiece], 'w')
    end
    startco = 1
    endco = seq.length
    index[:current].puts ">#{name}_#{startco}_#{endco}\n#{seq}"
    index[:length] += seq.size
    index[:index][name] = index[:num]
    index[:output].puts "#{name}\t#{startco}\t#{endco}\t#{index[:layer]}\t#{File.basename(index[:refpiece])}"
    
  else	# need to break up the ref seq
    0.upto((seq.length.to_f/index[:cut]).to_i) do |i|
      startco = i * index[:cut] + 1
      endco = startco + index[:cut] + index[:overlap]
      
      if endco > seq.length
        endco = seq.length
      end
      
      index[:current].close	if index[:current] != nil
      index[:num] += 1
      index[:length] = 0
      if index[:num] % 1000 == 1  # 1000 per dir
        index[:layer] +=1 
        laydir = index[:basedir] + "#{index[:layer]}" + '/'
        system("mkdir #{laydir}")
      end
      
      index[:refpiece] = index[:basedir] + "#{index[:layer]}" + '/' + 'piece_' + "#{index[:num]}" + '.fa'
      
      $stderr.puts index[:refpiece]
      
      index[:current] = File.new(index[:refpiece], 'w')
      ss = seq.slice(startco - 1, endco - startco + 1)
      tmpName = "#{name}_#{startco}_#{endco}"
      index[:current].puts ">#{tmpName}\n#{ss}"
      index[:length] += seq.size
      index[:index][tmpName] = index[:num]
      index[:output].puts "#{name}\t#{startco}\t#{endco}\t#{index[:layer]}\t#{File.basename(index[:refpiece])}"
    end
  end	
end

name = ''
seq = ''


File.new(ref, 'r').each do |line|
  if line=~ /^>(\S+)/
    formatDB(name, seq, $index)
    name = $1
    seq = ''
    #    $stderr.puts name
  else
    seq << line.chomp
  end
end
formatDB(name, seq, $index)
$index[:output].close
