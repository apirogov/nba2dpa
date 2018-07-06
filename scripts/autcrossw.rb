#!/usr/bin/env ruby
#transparent autcross wrapper that generates csv on the fly (appends after each input automaton)
#and also uses /dev/shm (RAMfs) for temp files
require 'open3'

SCRIPTPATH=`echo -n $( cd "$(dirname "$0")" ; pwd -P )`
CWD=`realpath .`.strip
SEED=`date +%s`.strip
DATE=`date +%Y%m%d-%H%M%S`.strip

#autcross with optional global timeout (the sanity checks can take long too)
tout = ARGV[0].to_i
ARGV.shift if tout > 0

tmpcsv=`mktemp -p /dev/shm`.strip
tmphoa=`mktemp -p /dev/shm`.strip
Signal.trap("INT") { `rm #{tmpcsv} #{tmphoa}`; exit }
Signal.trap("TERM") { `rm #{tmpcsv} #{tmphoa}`; exit }

csvfile=''
hoafile=''
ARGV.map! do |arg|
  if arg =~ /^--csv=/
    csvfile = arg.split('=')[1].strip
    arg = "--csv=#{tmpcsv}"
  end
  if arg =~ /^--save-bogus=/
    hoafile = arg.split('=')[1].strip
    arg = "--save-bogus=#{tmphoa}"
  end
  arg
end

cmd = "autcross '#{ARGV[0..-1].to_a.join("' '")}'"
cmd = "timeout #{tout} "+cmd if tout>0

# puts cmd
# exit

def each_aut(file, &fun)
  aut = ""
  file.each_line do |l|
    aut += l
    l.strip!
    aut = "" if l == "--ABORT--"
    if l == "--END--"
      fun.call(aut)
      aut = ""
    end
  end
end

open(csvfile, 'w') { |f| f.puts DATA.read } if csvfile!=''
each_aut(STDIN) do |a|
  Open3.popen3({"SPOT_TMPDIR" => '/dev/shm'}, cmd) do |stdin, stdout, stderr, wait_thr|
    stdin.write a
    stdin.close

    Thread.new do
      stdout.each_line do |l|
        puts l
      end
    end
    Thread.new do
      stderr.each_line do |l|
        STDERR.puts l
      end
    end

    wait_thr.join
    # exit_status = wait_thr.value.exitstatus
  end
  open(csvfile, 'a'){ |f| f.write `tail -n +2 #{tmpcsv}` }
  open(hoafile, 'a'){ |f| f.write open(tmphoa, 'r').read() }
end

`rm #{tmpcsv} #{tmphoa}`

__END__
"input.source","input.name","input.ap","input.states","input.edges","input.transitions","input.acc_sets","input.scc","input.nondetstates","input.nondeterministic","input.alternating","tool","exit_status","exit_code","time","output.ap","output.states","output.edges","output.transitions","output.acc_sets","output.scc","output.nondetstates","output.nondeterministic","output.alternating"
