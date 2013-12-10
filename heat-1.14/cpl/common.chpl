
const FPOPS_PER_POINT_PER_ITERATION = 12;

record params {
    var N : int;
    var M : int;
    var maxiter : int;
    var period : int;
    var threshold : real;
    var io_tmin : real;
    var io_tmax : real;
    var nthreads : int;
    var tinit;
    var tcond;
}

record results {
    var niter : int;
    var tmin : real;
    var tmax : real;
    var maxdiff : real;
    var tavg : real;
    var time : real;
}

proc print_header(){ 
    writeln("Output from heat-1.15 (pcs2013@list.uva.nl):\n\n",
               "Iterations",
               "T(min)     ", "T(max)     ", "T(diff)    ", "T(avg)     ", "Time      ", "FLOP/s    ");
}

proc report_results(p : params, r : results) {

    writeln(r.niter, " ",
            r.tmin, " ",
            r.tmax, " ",
            r.maxdiff, " ",
            r.tavg, " ",
            r.time, " ",
            p.N:real * p.M:real * 
            (r.niter:real * FPOPS_PER_POINT_PER_ITERATION +
                    r.niter:real / p.period:real) / r.time);
}
