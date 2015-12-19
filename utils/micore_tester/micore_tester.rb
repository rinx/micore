#!/usr/bin/env ruby
# -*- coding: utf-8 -*-

#
# micore_tester.rb
#
# Author: Rintaro Okamura
#
# Description:
#   This is a tester script for MICO-RE.
#   It makes a 2d histogram plot with differences from the true tau and cder.
#
#   RTM_BIN: Radiative transfer model should be prepared.
#   It should be given two parameters (tau and cder) as arguments
#   and it should provide reflectances for 2 wavelengths as stdout.
#     $ #{RTM_BIN} #{tau} #{cder}
#       REF1: #{reflectance1}
#       REF2: #{reflectance2}
#
# Dependency:
#   NArray gem
#   gnuplot gem
#
# Changelog:
#   20151219: First version
#

require 'narray'
require 'gnuplot'
require 'open3'
require 'optparse'

MICORE_BIN = '../../src/micore'
RTM_BIN    = 'sh ./rstar_wrapper.sh'

TAU_MIN  = 0.3
TAU_MAX  = 100.0
CDER_MIN = 5.0
CDER_MAX = 32.0

DEF_N_SAMPLE = 100
BIN_RESOLUTION = 20

params = ARGV.getopts('o:', 'p:', 'n:')

if ARGV.size < 1 then
  STDERR.puts <<heredoc
Usage: ruby #{$0} [-o txtname] [-p picname] [-n number] lutfile
  [-o txtname]: output txt filename.
  [-p picname]: output png filename (png, pdf, or tex).
  [-n number]: # of sample (default: #{DEF_N_SAMPLE}).
  lutfile: LUT binary file.
heredoc
  exit
end

# https://github.com/rdp/ruby_gnuplot/blob/master/lib/gnuplot.rb#L358
class SPlot < Gnuplot::Plot

  def initialize (io = nil, cmd = "splot")
    @cmd = cmd
    @settings = []
    @arbitrary_lines = []
    @data = []
    @styles = []
    yield self if block_given?
    $stderr.puts "writing this to gnuplot:\n" + to_gsplot + "\n" if $VERBOSE

    if io
      io << to_gsplot
      io << store_datasets
    end
  end

  def to_gsplot (io = "")
    @settings.each do |setting|
      io << setting.map(&:to_s).join(" ") << "\n"
    end
    @styles.each{|s| io << s.to_s << "\n"}
    @arbitrary_lines.each{|line| io << line << "\n" }

    io
  end

  def store_datasets (io = "")
    if @data.size > 0
      io << @cmd << " " << @data.collect { |e| e.plot_args }.join(", ")
      io << "\n"

      v = @data.collect { |ds| ds.to_gsplot }
      io << v.compact.join("e\n")
    end
  end
end

class Array
  def to_gsplot
    f = ""

    if ( self[0].kind_of? Array ) then
      x = self[0]
      y = self[1]
      d = self[2]

      x.each_with_index do |xv, i|
        y.each_with_index do |yv, j|
          f << [ xv, yv, d[i][j] ].join(" ") << "\n"
        end
        f << "\n"
      end
    elsif ( self[0].kind_of? Numeric ) then
      self.length.times do |i| f << "#{self[i]}\n" end
    else
      self[0].zip( *self[1..-1] ).to_gsplot
    end

    f
  end
end

lutfilepath = ARGV[0]
outfilepath = params['o']
picfilepath = params['p']
picfilepath ||= "none"
nsample = params['n']
nsample ||= DEF_N_SAMPLE
nsample = nsample.to_i

# check LUT file
unless File.exist?(File.expand_path(lutfilepath)) then
  STDERR.puts "LUT file #{lutfilepath} doesn't exist."
  exit
end

# initialize
obs_tau  = NArray.sfloat(nsample)
obs_cder = NArray.sfloat(nsample)
ref1     = NArray.sfloat(nsample)
ref2     = NArray.sfloat(nsample)
est_tau  = NArray.sfloat(nsample)
est_cder = NArray.sfloat(nsample)
est_cost = NArray.sfloat(nsample)

# generate random cloud properties
STDOUT.puts "Generating random taus and cders..."
obs_tau.random!(TAU_MAX - TAU_MIN)
obs_cder.random!(CDER_MAX - CDER_MIN)
obs_tau  = obs_tau + TAU_MIN
obs_cder = obs_cder + CDER_MIN

nsample.times do |i|
  STDOUT.puts "RTM calculation #{i}: started"
  # RTM calculation
  rtm_stdout, rtm_status = Open3.capture2("#{RTM_BIN} #{obs_tau[i]} #{obs_cder[i]}")
  STDOUT.puts "RTM calculation #{i}: finished"

  # parse RTM output
  rtm_stdout.each_line do |line|
    case line
    when /^\s*REF1:/
      ref1[i] = line.gsub(/^\s*REF1:\s*/,'').to_f
    when /^\s*REF2:/
      ref2[i] = line.gsub(/^\s*REF2:\s*/,'').to_f
    end
  end

  # MICO-RE calculation
  STDOUT.puts "MICO-RE calculation #{i}: started"
  micore_stdout, micore_status = Open3.capture2("#{MICORE_BIN} #{lutfilepath} #{ref1[i]} #{ref2[i]}")
  STDOUT.puts "MICO-RE calculation #{i}: finished"

  # parse MICO-RE output
  micore_stdout.each_line do |line|
    case line
    when /^\s*TAU:/
      est_tau[i]  = line.gsub(/^\s*TAU:\s*/,'').to_f
    when /^\s*CDER:/
      est_cder[i] = line.gsub(/^\s*CDER:\s*/,'').to_f
    when /^\s*COST:/
      est_cost[i] = line.gsub(/^\s*COST:\s*/,'').to_f
    end
  end
end

# calculate differences between estimation and observation
STDOUT.puts "Calculating 2d histogram"
diff_tau  = est_tau  - obs_tau
diff_cder = est_cder - obs_cder

absmax_tau  = diff_tau.abs.max
absmax_cder = diff_cder.abs.max

# calculate 2d histogram
nx = BIN_RESOLUTION + 1
ny = BIN_RESOLUTION + 1
interval_x = (2 * absmax_tau  / BIN_RESOLUTION)
interval_y = (2 * absmax_cder / BIN_RESOLUTION)
dx = interval_x / 2
dy = interval_y / 2
x = []
y = []
(BIN_RESOLUTION + 1).times do |i|
  x.push(absmax_tau  * (-1) + interval_x * i)
  y.push(absmax_cder * (-1) + interval_y * i)
end

hist = Array.new(BIN_RESOLUTION + 1).map{Array.new(BIN_RESOLUTION + 1, 0)}
nx.times do |ix|
  ny.times do |iy|
    insidebin_tau_index  = ((diff_tau  >= (x[ix] - dx)).where.to_a & (diff_tau  < (x[ix] + dx)).where.to_a)
    insidebin_cder_index = ((diff_cder >= (y[iy] - dy)).where.to_a & (diff_cder < (y[iy] + dy)).where.to_a)
    hist[ix][iy] = (insidebin_tau_index & insidebin_cder_index).length
  end
end

# write output
if outfilepath
  File.open(outfilepath, "w") do |f|
    f.puts "# obs_tau, obs_cder, ref1, ref2, est_tau, est_cder, est_cost, diff_tau, diff_cder"
    nsample.times do |i|
      f.puts "#{obs_tau[i]}, #{obs_cder[i]}, #{ref1[i]}, #{ref2[i]}, #{est_tau[i]}, #{est_cder[i]}, #{est_cost[i]}, #{diff_tau[i]}, #{diff_cder[i]}"
    end
  end
end

# make plot
Gnuplot.open do |gp|
  SPlot.new(gp) do |plot|
    STDOUT.puts "plotting..."

    case picfilepath
    when /^none$/
      plot.terminal 'x11'
    when /\.tex$/
      plot.terminal "tikz"
      plot.output File.expand_path(picfilepath)
    when /\.pdf$/
      plot.terminal "pdf"
      plot.output File.expand_path(picfilepath)
    when /\.png$/
      plot.terminal "png"
      plot.output File.expand_path(picfilepath)
    else
      STDERR.puts "Error: unknown filetype"
      exit
    end

    plot.title "differences between estimation and observation"

    plot.xlabel "tau(est) - tau(obs)"
    plot.ylabel "cder(est) - cder(obs)"
    plot.cblabel "frequency"

    plot.set "key off"

    plot.set "pm3d map"
    plot.set "size square"

    plot.xrange "[-#{absmax_tau}:#{absmax_tau}]"
    plot.yrange "[-#{absmax_cder}:#{absmax_cder}]"

    plot.cbrange "[0:#{hist.flatten.max}]"

    plot.set "palette defined ( 0 '#000090', 1 '#000fff', 2 '#0090ff', 3 '#0fffee', 4 '#90ff70', 5 '#ffee00', 6 '#ff7000', 7 '#ee0000', 8 '#7f0000')"

    plot.data <<
    Gnuplot::DataSet.new([x, y, hist]) do |ds|
    end

    STDOUT.puts "finish to plot."
  end
end


