## given a SNP list from the whole genome and targeted genomic region, grab all the SNPs in the targeted region

snplist = ARGV[0]
targeted = ARGV[1]

$gray = 200
$ranges = {}

File.new(targeted, 'r').each do |line|
  if line=~ /^(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)/
    tag, chr, offset, s,e =$1, $2, $3.to_i, $4.to_i, $5.to_i
    ref = "chr" + chr
    $ranges[ref] = [] unless $ranges.key?(ref)

    $ranges[ref] << Range.new(s - $gray, e + $gray)
  end
end

# $ranges.each_key do |ref|
#  $stderr.puts "#{ref} \t #{$ranges[ref].join("\t")}"
# end

File.new(snplist,'r').each do |line|
  if line=~ /^(chr\S+)\s+(\d+)/
    ref , pos = $1, $2.to_i
    if $ranges.key?(ref)
#      $stderr.puts ref, pos
      $ranges[ref].each do |r|
        if r.include?(pos) # the SNP lies in the targeted region

          puts line.chomp
        end
      end
    end
  end
end
