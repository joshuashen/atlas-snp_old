## given a SNP list from the whole genome and targeted genomic region, grab all the SNPs in the targeted region 

require 'getoptlong'

opts = GetoptLong.new(
    ["--input", "-i", GetoptLong::REQUIRED_ARGUMENT],
    ["--target", "-t", GetoptLong::REQUIRED_ARGUMENT],
    ["--output", "-o", GetoptLong::OPTIONAL_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

if optHash.key?("--help") or !optHash.key?("--input") or !optHash.key?("--target")
  $stderr.puts "Usage: ruby __.rb -i input_SNP_list -t target_regions [ -o prefix_of_output ]"
  exit
end

if optHash.key?("--output")
  prefix = optHash["--output"]
else
  prefix = optHash["--target"]
end

$ranges = {} # store the genomic locations of all genes
$snps = {}  # store the number of SNPs for each gene ddb

# ddb name (as key) -> gene name (as value),  a gene could have multiple ddb names; a ddb always points to one gene.
$ddbh = {}

$chromosome_change={"M"=>"DDB0169550","2F"=>"DDB0215018", "BF"=>"DDB0220052", "3F"=>"DDB0215151", "1"=>"DDB0232428", "2"=>"DDB0232429","3"=>"DDB0232430", "4"=>"DDB0232431", "5"=>"DDB0232432", "6"=>"DDB0232433"}


snpout = File.new(prefix + ".SNP_list", "w")
statout = File.new(prefix + ".SNP_stats", "w")

# read the file containing gene information
File.new(optHash["--target"], 'r').each do |line|
  if line=~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(.*)$/
    gene, ddb, chr, s, e , strand, desc = $1,$2,$3, $4.to_i, $5.to_i, $6, $7
    ref = $chromosome_change[chr]

    # initiate the hash for this ref if not already
    $ranges[ref] = Hash.new unless $ranges.key?(ref)

    # a hash for this ddb, ddb name as the key, Range as the value
    $ranges[ref][ddb] = Range.new(s,e)

    $ddbh[ddb] = gene
    $snps[ddb] = 0
  end
end

# read SNP file                                                                                       
File.new(optHash["--input"],'r').each do |line|
  if line=~ /^(\S+)\s+(\d+)/
    ref , pos = $1, $2.to_i
    
    # an array to store the ddb's of all relevant genes                                               
    genes = []
    if $ranges.key?(ref)
      # go through all genes (by ddb) in this chromosome                                              
      $ranges[ref].each_key do |ddb|
        # judge if this number "pos" is within the Range $ranges[ref][ddb]                            
        if $ranges[ref][ddb].include?(pos)
          genes << ddb
          $snps[ddb] += 1
        end
      end
    end
    
    #get at least one gene                                                                            
    if genes.size > 0
      gnames = {}
      genes.each do |ddb|
        gnames[$ddbh[ddb]] = 1
      end

      snpout.puts gnames.keys.join(",") + "\t" + genes.join(",") + "\t" + line.chomp!
      
    end
  end
end

# print the number of SNPs for each gene, a ddb a line                                                
$snps.each_key do |ddb|
  statout.puts "#{ddb}\t#{$ddbh[ddb]}\t#{$snps[ddb]}"
end

snpout.close
statout.close
