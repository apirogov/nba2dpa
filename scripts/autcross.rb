#!/usr/bin/env ruby
require 'open3'

SCRIPTPATH=`echo -n $( cd "$(dirname "$0")" ; pwd -P )`
CWD=`realpath .`.strip
SEED=`date +%s`.strip
DATE=`date +%Y%m%d-%H%M%S`.strip

def wrap(str); 'cat %H | ' + str + ' > %O'; end

def oneof(*a); a.map(&->(x){[x]}); end

def prefs(*a); 0.upto(a.length).collect(&->(i){a[0...i]}); end

#list monad
class Array; def prod(rhs); self.product(rhs).map!(&:flatten); end; end

def translator_args(tarr); tarr.map(&->(s){"-t '#{wrap s}'"}).join(' '); end

def compile(a); a.map(&->(l){l.map(&->(x){"-"+x}).join(" ")}); end
#e.g.:
#compile oneof("a","b").prod prefs("c","d","e")

# ----
# compact representation of different argument combinations
nbadet_argsets = compile([['k','j','t','n','a','d']].prod(
                            oneof("u0","u1","u2").prod(
                              prefs("e","b","l","m","c")
                            )
                          )
                        )

translators = ["autfilt -D -C -P"]
nbadet_argsets.each do |s|
  translators.push("#{SCRIPTPATH}/../build/bin/nbadet #{s}")
end

output_args = "--save-bogus=failed_#{DATE}.hoa --csv=stats_#{DATE}.csv"
cmd = "./autcrossw.rb #{translator_args translators} #{output_args} #{ARGV.join(' ')}"
# puts cmd
# exit

# ----
# call autcross, pass through IO
Open3.popen3(cmd) do |stdin, stdout, stderr, wait_thr|
  # consume outputs in STDOUT and STDERR, otherwise the process
  # may be blocked if it produces a lot of outputs
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

  while data = STDIN.read(256) # read output of the upstream command
    stdin.write(data)          # manually pipe it to the ffmpeg command
  end
  stdin.close
  wait_thr.join
end
