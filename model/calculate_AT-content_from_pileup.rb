## input: pileup, window size

# ignore the very edges of chromosomes

$wsize = 200
$bases = {'A' =>0, 'C' => 0, 'G' => 0, 'T' => 0, 'N' => 0}
$head = nil
$window = []


def main
  p, $wsize = getArgs()
  middle = $wsize / 2 
#  $stderr.puts middle

  File.new(p,'r').each do |line|
    cols = line.split(/\s+/)
    base = {}
    base[:pos], base[:rbase], base[:cov] = cols[1], cols[2], cols[3]
    base[:at] = 0.5
    $window << base
    $bases[base[:rbase]] += 1

    if $window.size > $wsize  # compute the AT of the middle base
      base = $window[middle]
      $head = $window.shift
      $bases[$head[:rbase]] -= 1
      base[:at] = (100 * ( $bases['A'] + $bases['T'] ) / $wsize.to_f).round/100.0 # update AT content, only keep two decimals
      puts "#{$head[:pos]}\t#{$head[:rbase]}\t#{$head[:cov]}\t#{$head[:at]}"
    end
  end
end

def getArgs
  pileup = ARGV[0]
  wsize = ARGV[1]
  
  if wsize == nil
    wsize = 230
  else
    wsize = wsize.to_i
  end
  return [pileup, wsize]
end


main()
