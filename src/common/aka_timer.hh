/**
 * @file   aka_timer.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Tue Jun 17 2014
 *
 * @brief  timer object
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TIMER_HH__
#define __AKANTU_TIMER_HH__


#include <cmath>
#include <limits>
#include <ctime>
#include <sys/time.h>


__BEGIN_AKANTU__


//  time_t       tv_sec;   /* seconds since Jan. 1, 1970 */
//  suseconds_t  tv_usec;  /* and microseconds */

//  time_t -> long
//  suseconds_t -> __int32_t



using std::cout;
using std::endl;

class time {
  
  typedef std::time_t time_t;
  
  time_t microseconds_;
  time_t seconds_;
  time_t minutes_;
  time_t hours_;
  time_t days_;
  
public:
  

  explicit time(const timeval& t)
  : microseconds_(t.tv_usec), seconds_(t.tv_sec), minutes_(), hours_(), days_() {
    adjust();
  }

  
  explicit time(time_t ms = 0, time_t s = 0, time_t min = 0, time_t h = 0, time_t d = 0)
  : microseconds_(ms), seconds_(s), minutes_(min), hours_(h), days_(d) {
    adjust();
  }
  
  time& operator+=(const timeval& t) {
    microseconds_ += t.tv_usec;
    seconds_ += t.tv_sec;
    adjust();
    return *this;
  }
  
  time& operator+=(const time& t) {
    microseconds_ += t.microseconds_;
    seconds_ += t.seconds_;
    minutes_ += t.minutes_;
    hours_ += t.hours_;
    days_ += t.days_;
    adjust();
    return *this;
  }
  
  time operator+(const time& t2) {
    time t(*this);
    t += t2;
    return t;
  }
  
  friend std::ostream& operator<<(std::ostream& os, const time& t) {
    if (t.days_ != 0)
      os<<t.days_<<" d, ";
    if (t.hours_ != 0)
      os<<t.hours_<<" h, ";
    if (t.minutes_ != 0)
      os<<t.minutes_<<" m, ";
    if (t.seconds_ != 0)
      os<<t.seconds_<<" s, ";
    os<<t.microseconds_<<" µs";
    return os;
  }
  
private:
  
  void adjust() {

    if (microseconds_ >= 1000000) {
      seconds_ += microseconds_ / 1000000;
      microseconds_ = microseconds_ % 1000000;
    }
    if (seconds_ >= 60) {
      minutes_ += seconds_ / 60;
      seconds_ = seconds_ % 60;
    }
    if (minutes_ >= 60) {
      hours_ += minutes_ / 60;
      minutes_ = minutes_ % 60;
    }
    if (hours_ >= 24) {
      days_ += hours_ / 24;
      hours_ = hours_ % 24;
    }
  }
};


///////////////////////////////////////////////////////////////////////////////
// timer classes

template <typename timeval>
void subtract(const timeval& x, timeval y, timeval& r) {
    
  // perform the carry for the later subtraction by updating y
  if (x.tv_usec < y.tv_usec) {
    int nsec = (y.tv_usec - x.tv_usec) / 1000000 + 1;
    y.tv_usec -= 1000000 * nsec;
    y.tv_sec += nsec;
  }
  if (x.tv_usec - y.tv_usec > 1000000) {
    int nsec = (x.tv_usec - y.tv_usec) / 1000000;
    y.tv_usec += 1000000 * nsec;
    y.tv_sec -= nsec;
  }
  
  // compute the time remaining to wait
  // tv_usec is certainly positive
  r.tv_sec = x.tv_sec - y.tv_sec;
  r.tv_usec = x.tv_usec - y.tv_usec;
}


class ctimer {
  
  timeval t1_;
  
public:

  ctimer()
  { reset(); }
    
  timeval tac() const {
    timeval t2;
    gettimeofday(&t2, NULL);
    timeval r;
    subtract(t2, t1_, r);
    return r;
  }
  
  inline void reset()
  { gettimeofday(&t1_, NULL); }
  
  friend std::ostream& operator<<(std::ostream& os, const ctimer& t) {
    cout<<time(t.tac());
    return os;
  }
};


class cpu_timer {
  
  std::clock_t t1_;
  
public:
  cpu_timer() : t1_(std::clock()) {}
  
  
  std::clock_t tac() const {
    return std::clock() - t1_;
  }
  
  inline void reset() { t1_ = std::clock(); }
  
  inline double seconds()
  { return static_cast<double>(tac())/CLOCKS_PER_SEC; }
  
  friend std::ostream& operator<<(std::ostream& os, const cpu_timer& t) {
    os<<t.tac();
    return os;
  }
};


__END_AKANTU__

#endif /* __AKANTU_TIMER_HH__ */