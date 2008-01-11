## 
blat = ARGV[0]
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


File.new(blat,'r').each do |line|
  cols = line.split(/\s+/)
  ref , s ,e  = cols[2], cols[7].to_i, cols[8].to_i
  if $ranges.key?(ref)
#      $stderr.puts ref, pos                                                                          
    $ranges[ref].each do |r|
      if (s >= r.first and s <= r.last ) or ( e <= r.last and e >= r.first)
        puts line.chomp
      end
    end
  end
end

