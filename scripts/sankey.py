#!/usr/bin/python
r""" 
    Plot the Orion proto transition matrix for a certain level
"""

import os
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from operator import itemgetter
from argparse import RawTextHelpFormatter

def printChildren(pref, level):
    for iP in xrange(nprotos):
        label = "%s %d.%d" % (pref, iP, level)
        if label in paths and paths[label] >= minfreq:
            print "      ['%s', '%s', %d]," % (pref, label, paths[label])
            printChildren(label, level+1)

def printParents(pref, level):
    for iP in xrange(nprotos):
        label = "%s %d.%d" % (pref, iP, level)
        if label in paths and paths[label] >= minfreq:
            print "      ['%s', '%s', %d]," % (label, pref, paths[label])
            printParents(label, level+1)

def printHeader():
    print """<html>
<body>
<script type="text/javascript"
   src="https://www.google.com/jsapi?autoload={'modules':[{'name':'visualization','version':'1.1','packages':['corechart', 'sankey']}]}">
</script>
 
<script type="text/javascript">
  google.setOnLoadCallback(drawChart);

  function drawChart() {
"""

def printGlobalHeader():
    print """<html>
<body>
<script type="text/javascript"
   src="https://www.google.com/jsapi?autoload={'modules':[{'name':'visualization','version':'1.1','packages':['sankey']}]}">
</script>

<div id="sankey_multiple" style="width: 1400px; height: 800px;"></div>

<script type="text/javascript">
  google.setOnLoadCallback(drawChart);
    function drawChart() {
      var data = new google.visualization.DataTable();
      data.addColumn('string', 'From');
      data.addColumn('string', 'To');
      data.addColumn('number', 'Weight');
      data.addRows([
""",
    
def printGlobalFooter():
    print """
    ]);

    // Set chart options
    var options = {
      width: 1300,
      sankey: { iterations: 128},
    };

    // Instantiate and draw our chart, passing in some options.
    var chart = new google.visualization.Sankey(document.getElementById('sankey_multiple'));
    chart.draw(data, options);
  }
</script>
</body>
"""

def get_args():
    r""" Parse arguments for the program    
    """
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument("pathsfile", help="File containing the proto paths of each sequence.")
    parser.add_argument("ftrsfile", nargs='?', default=None, help="File containing the key features for each proto (required for modes proto and rproto).")
    parser.add_argument("-f","--minfreq", type=int, default=1, help="Minimum frequency to include a transition in the chart.")
    parser.add_argument("-l","--maxlevel", type=int, default=4, help="Maximum number of levels to show in the chart (for mode global).")
    parser.add_argument("-m", "--mode", type=str, default="global", help="""Type of chart to show:
    global  Show global evolution chart [default]
    proto   Show per-proto evolution
    rproto  Show reverse per-proto evolution
"""
)
    
    return parser, parser.parse_args()


if __name__ == '__main__':

    parser, args = get_args()  # command line arguments
    
    pathsfile = args.pathsfile
    ftrsfile = args.ftrsfile
    minfreq = args.minfreq
    nlevels = 0
    nprotos = 0
    nseq = 0
    
    if not os.path.isfile(pathsfile):
        print "Missing required file: pathsfile. '%s' is not a file.\n" % pathsfile
        parser.print_help()
        sys.exit()
        
    if args.mode != "global" and (ftrsfile is None or not os.path.isfile(ftrsfile)):
        print "Missing required file for mode %s: ftrsfile. '%s' is not a file.\n" % (args.mode, ftrsfile)
        parser.print_help()
        sys.exit()
        
    # read in the paths
    paths = {}
    counts = {}
    with open(pathsfile, "r") as fh:
        for line in fh:
            if not line:
                continue
            fields = line.strip().split(" ")
            prev = 'S'
            if len(fields) > nlevels:
                nlevels = len(fields)
                
            if args.mode == "rproto":
                fields.reverse()
                
            for i,f in enumerate(fields):
                f = int(f)
                if f > nprotos:
                    nprotos = f
                    
                if prev not in counts:
                    counts[prev] = 0
                counts[prev] += 1
                    
                if args.mode == "global":
                    label = "%s P%d-%d" % (prev, f, i+1)
                    prev = "P%d-%d" % (f, i+1)
                else:
                    if i == 0:
                        label = "%d.%d" % (f, i)
                    else:
                        label = "%s %d.%d" % (prev, f, i)
                    prev = label
                if label not in paths:
                    paths[label] = 0
                paths[label] += 1
                
            if prev not in counts:
                counts[prev] = 0
            counts[prev] += 1
            
            nseq += 1
    
    nprotos += 1
    
    # Read the features 
    vcounts = {}
    norms = {}
    sqes = {}
    nftrs = {}
    pwgts = {}
    cwgts = {}
    labels = {}
    if args.mode == "proto" or args.mode == "rproto":
        with open(ftrsfile, "r") as fh:
            for line in fh:
                if not line:
                    continue
                dt = line.strip().split("\t")
                if len(dt) < 2:
                    raise ValueError("Invalid data in line. Metadata should be separated by a tab from the data.\n %s" % line)
                fields = dt[0].split()
                if len(fields) < 4:
                    raise ValueError("Invalid metadata in line. Metadata should contain pid, pcount, pint, and psqe.\n %s" % line)
                pid = int(fields[0])
                vcounts[pid] = int(fields[1])
                norms[pid] = float(fields[2])
                sqes[pid] = float(fields[3])
                fields = dt[1].split()
                if len(fields) % 5 != 0:
                    raise ValueError("Invalid data in line. Data should be a set of (fint cint fdis sum flabel) tuples, separated by spaces.\n %s" % line)
                nftrs[pid] = len(fields)/5
                
                rp = 0.0
                rc = 0.0
                pwgts[pid] = []
                cwgts[pid] = []
                labels[pid] = []
                for i in xrange(nftrs[pid]):
                    pwgts[pid].append(float(fields[5*i+0]))
                    cwgts[pid].append(float(fields[5*i+1]))
                    labels[pid].append(fields[5*i+4])
                    rp += pwgts[pid][i]
                    rc += cwgts[pid][i]
                        
                # Add an entry for the 'rest'
                pwgts[pid].append( 1.0-rp )
                cwgts[pid].append( 1.0-rc )
                labels[pid].append("all other")
                nftrs[pid] += 1
                
                
    if args.mode == "global":
        
        printGlobalHeader()
        for label in paths.keys():
            pfrom, pto = label.split()
            _, level = pto.split('-')
            if int(level) > args.maxlevel:
                continue
            if paths[label] >= args.minfreq:
                print "          ['%s (%04d)', '%s (%04d)', %d]," % (pfrom, counts[pfrom], pto, counts[pto], paths[label])
        printGlobalFooter()
        
    elif args.mode == "proto" or args.mode == "rproto":
        
        printHeader()
        for iP in xrange(nprotos):
            pref = "%d.0" % iP
            if pref in paths and paths[pref] >= args.minfreq:
                print "    var data%d = new google.visualization.DataTable();" % (iP)
                print "    data%d.addColumn('string', 'From');" % (iP)
                print "    data%d.addColumn('string', 'To');" % (iP)
                print "    data%d.addColumn('number', 'Weight');" % (iP)
                print "    data%d.addRows([" % (iP)
                if args.mode == "proto":
                    print "      ['S%d', '%s', %d]," % (paths[pref], pref, paths[pref])
                    printChildren(pref, 1)
                else: 
                    print "      ['%s', 'E%d', %d]," % (pref, paths[pref], paths[pref])
                    printParents(pref, 1)
                print "    ]);"
            
                print "    var cdata%d = google.visualization.arrayToDataTable([" % (iP)
                print "      ['App', 'Usage'],"
                
                # features data
                for iF in xrange(nftrs[iP]):
                    print "      ['%s', %.2f]," % (labels[iP][iF], 100.0*cwgts[iP][iF])
                print "    ]);"
            
                print "    var pdata%d = google.visualization.arrayToDataTable([" % (iP)
                print "      ['App', 'Usage'],"
                
                # features data
                for iF in xrange(nftrs[iP]):
                    print "      ['%s', %.2f]," % (labels[iP][iF], 100.0*pwgts[iP][iF])
                print "    ]);"
                
                print "    var soptions%d = { height: 400, sankey: { iterations: 128}, };" % (iP)
                print "    var chart%d = new google.visualization.Sankey(document.getElementById('sankey%d'));" % (iP, iP)
                print "    chart%d.draw(data%d, soptions%d);\n" % (iP, iP, iP)
            
                print "    var diffChart%d = new google.visualization.PieChart(document.getElementById('proto%d'));" % (iP, iP)
                print "    var diffData%d = diffChart%d.computeDiff(cdata%d, pdata%d);" % (iP, iP, iP, iP)
                print "    var poptions%d = { pieSliceText: 'none' };" % (iP)
                print "    diffChart%d.draw(diffData%d, poptions%d);\n" % (iP, iP, iP)
        print """ }\n</script>\n\n"""
        
        for iP in xrange(nprotos):
            print "<table>"
            print "  <tr><td>"
            print "    <table>"
            print "      <tr><td>Proto %d: #Weeks: %d, AvgUsage: %.1f, Dist2Center: %.1f</td></tr>" % (iP, vcounts[iP], norms[iP], sqes[iP])
            print "      <tr><td><div id='proto%d' style='width:600px; height: 430px; display: inline-block;'></div></td></tr>" % (iP)
            print "    </table>"
            print "  </td><td>"
            print "    <div id='sankey%d' style='width: 900px; height: 450px;'></div>" % (iP)
            print "  </td></tr>"
            print "</table><hr>\n"
        
        print """</body>"""
        
    else:
        
        print "Invalid mode."
        parser.print_help()
        sys.exit()
        