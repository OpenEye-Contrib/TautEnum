// chrono.h - A Chronograph class
//            for program timing.
// Greg Messer

/*
 NOTES:

 This class provides a Chronograph for an
 application.  It works like a stopwatch.
 It provides timing for laps (time since
 last lap) and splits (same as elapsed time).
 It can be stopped, restarted and reset to
 zero.

 This class uses the ANSI C clock() function
 to get the time.  The clock() function
 returns a clock_t that is the number of
 clock ticks since a program started.
 The lap, split and elapsed times are
 calculated as the difference between two
 clock_t values.  If the program is a
 long-running one, the clock_t may overflow.
 So, don't use this for timing longer than
 the value of MAX_LONG divided by the ANSI C
 macro CLOCKS_PER_SECOND.

 The compiler-supplied copy constructor and
 assignment operator are OK for this class,
 so they are not implemented here.

 Use the istream function precision(2) to
 format the returned double values from the
 member functions.

 Use "%.2lf" as the printf() format string
 to format the returned double values from
 the member functions.

 */

#ifndef CHRONO_H
#define CHRONO_H

#include <time.h>

class Chronograph
{
private:

   // when the Chronograph was started
   clock_t _start;

   // when the Chronograph was stopped
   // If _stopped != 0, the Chronograph
   // is stopped.
   clock_t _stopped;

   // when the last lap was read
   clock_t _lastlap;

   // the last elapsed time taken
   double _elapsed;

   /*
    diff() calculates the difference IN
    SECONDS between two clock_t values.
    */

   double diff(clock_t st, clock_t end)
   {
      return double(end - st) / CLOCKS_PER_SEC;
   }

public:

   // default constructor
   Chronograph()
   {
      reset();
      start();
   }

   // destructor
   ~Chronograph()
   {
   }

   /*
    start() either starts the Chronograph
    at the current time or restarts the
    Chronograph, adjusting the start and
    lastlap times by the amount of time that
    the chrono was stopped.
    Calling start() while the Chronograph is
    running restarts it at 0.0.
    */

   double start(void)
   {
      // current clock
      clock_t ct;

      // how long the Chronograph was stopped
      clock_t ctd;

      ct = clock();

      if(ct == 0)
         ct++;

      if(_stopped)
      {
         ctd = ct - _stopped;
         _start += ctd;
         _lastlap += ctd;

         _stopped = 0;

         return elapsed();
      }
      else
      {
         _lastlap = ct;
         _start = ct;
         _elapsed = 0.0;

         return _elapsed;
      }
   }

   /*
    stop() records the clock value when the
    Chronograph is stopped.  When stopped,
    the recorded value is used in the elapsed
    time calculation.  The Chronograph
    remains stopped until the start()
    function is called to restart it.
    */

   double stop(void)
   {
      if(!_stopped)
      {
         _stopped = clock();

         if(_stopped == 0)
            _stopped++;
      }

      return elapsed();
   }

   /*
    reset() initializes all times to the
    current time, and stops the Chronograph
    The start() function must be called to
    restart it.
    */

   double reset(void)
   {
      clock_t ct;     // current clock

      ct = clock();

      if(ct == 0)
         ct++;

      _lastlap = ct;
      _start = ct;
      _stopped = ct;
      _elapsed = 0.0;

      return _elapsed;
   }

   /*
    lap() returns the number of seconds
    since lap() was last called.
    */

   double lap(void)
   {
      clock_t ct;     // current clock
      double dl;      // lap time

      if(_stopped)
         ct = _stopped;
      else
      {
         ct = clock();

         if(ct == 0)
            ct++;
      }

      dl = diff(_lastlap, ct);

      if(!_stopped)
            _lastlap = ct;

      return dl;
   }

   /*
    split() returns the current split time,
    which is the elapsed time.  Split is a
    common stopwatch function so it is
    included here.
    */

   double split(void)
   {
      return elapsed();
   }

   /*
    isstopped() returns the value of the
    private _stopped member variable.  This
    is non-zero if the Chronograph is
    stopped, and zero if it is running.
    */

   int isstopped(void)
   {
      if(_stopped)
         return 1;
      else
         return 0;
   }

   /*
    elapsed() returns the elapsed time in
    seconds that the Chronograph has been
    running.
    */

   double elapsed(void)
   {
      clock_t ct;     // current clock

      if(_stopped)
         ct = _stopped;
      else
      {
         ct = clock();

         if(ct == 0)
            ct++;
      }

      _elapsed = diff(_start, ct);

      return _elapsed;
   }

   /*
    elapsedHMS() returns the elapsed time in
    seconds that the Chronograph has been
    running.  In addition, it breaks out
    the elapsed time into hours, minutes,
    and seconds.  Those are returned in the
    parameters passed in.
    */

   double elapsedHMS(double &Hours,
                     double &Mins,
                     double &Secs);
};

#endif

// end of chrono.h

