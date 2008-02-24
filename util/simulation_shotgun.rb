# simulate shotgun reads from reference genomes

require 'getoptlong'

opts = GetoptLong.new(
	["--input", "-i", GetoptLong::REQUIRED_ARGUMENT],
	["--coverage", "-c", GetoptLong::REQUIRED_ARGUMENT],
	["--readSize", "-s", GetoptLong::OPTIONAL_ARGUMENT],
	["--fixedSize", "-f", GetoptLong::OPTIONAL_ARGUMENT],
	["--help", "-h", GetoptLong::NO_ARGUMENT]
)


optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

if optHash.key?("--help") or (!optHash.key?("--readSize") and !optHash.key?("--fixedSize") )
	$stderr.puts "Usage: ruby #{$0} -i input.fa -c depth_coverage [-s read_size_distribution] [-f fixed_size] > simulated.fa"
	$stderr.puts "     -i  input.fa  -- the reference sequence targeted for shotgun sequencing"
	$stderr.puts "     -c  depth coverage"
	$stderr.puts "     -s  a file contain the list of real read size distribtution"
	$stderr.puts "     -f  fixed read size. Note: -s or -f must be provided; if -s is present, -f is ignored"
	exit
end	


$seq = {}
$sizes = []
$mean = 0
ref = ''
File.new(optHash["--input"], 'r').each do |line|
	if line=~ /^>(\S+)/
		ref = $1
		$seq[ref] = ''
	else
		$seq[ref] << line.chomp
	end
end 


if optHash.key?("--readSize")

	File.new(optHash["--readSize"], "r").each do |line|
		if line=~ /(\d+)/
			$sizes << $1.to_i
		end
	end	
	t = 0
	$sizes.map {|i| t+= i}
	$mean = t/$sizes.size.to_f
else
	sizes << optHash["--fixedSize"].to_i
	$mean = sizes[0]
end

sizesL = $sizes.size		

$seq.each_key do |ref|
	ss = $seq[ref].length  
	num = (ss * optHash["--coverage"].to_f / $mean).to_i
	starts = Array.new(num).map {rand(ss - 350)}
	ends = []
	starts.sort.each do |s|
		e = [s + $sizes[rand(sizesL)], ss].min
		puts ">#{ref}_#{s}_#{e}"
		puts $seq[ref][s-1..e-1]
	end
end

