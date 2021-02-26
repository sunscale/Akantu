/**
 * @file   aka_error.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jun 14 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  error management and internal exceptions
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include <set>
#include <sstream>
#include <typeinfo>
#include <utility>
#include <vector>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_ERROR_HH_
#define AKANTU_ERROR_HH_

namespace akantu {
/* -------------------------------------------------------------------------- */
enum DebugLevel {
  dbl0 = 0,
  dblError = 0,
  dblAssert = 0,
  dbl1 = 1,
  dblException = 1,
  dblCritical = 1,
  dbl2 = 2,
  dblMajor = 2,
  dbl3 = 3,
  dblCall = 3,
  dblSecondary = 3,
  dblHead = 3,
  dbl4 = 4,
  dblWarning = 4,
  dbl5 = 5,
  dblInfo = 5,
  dbl6 = 6,
  dblIn = 6,
  dblOut = 6,
  dbl7 = 7,
  dbl8 = 8,
  dblTrace = 8,
  dbl9 = 9,
  dblAccessory = 9,
  dbl10 = 10,
  dblDebug = 42,
  dbl100 = 100,
  dblDump = 100,
  dblTest = 1337
};

/* -------------------------------------------------------------------------- */
#define AKANTU_LOCATION                                                        \
  "(" << std::string(__func__) << "(): " << std::string(__FILE__) << ":"       \
      << std::to_string(__LINE__)                                              \
      << ")" // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)

/* -------------------------------------------------------------------------- */
namespace debug {
  void setDebugLevel(const DebugLevel & level);
  const DebugLevel & getDebugLevel();

  void initSignalHandler();
  std::string demangle(const char * symbol);
  template <class T> std::string demangle() {
    return demangle(typeid(T).name());
  }

  template <class T> std::string demangle(const T & t) {
    return demangle(typeid(t).name());
  }

  auto exec(const std::string & cmd) -> std::string;
  auto getBacktrace() -> std::vector<std::string>;
  void
  printBacktrace(const std::vector<std::string> & backtrace = getBacktrace());

  void exit(int status) __attribute__((noreturn));
  /* ------------------------------------------------------------------------ */
  /// exception class that can be thrown by akantu
  class Exception : public std::exception {
    /* ---------------------------------------------------------------------- */
    /* Constructors/Destructors                                               */
    /* ---------------------------------------------------------------------- */
  protected:
    explicit Exception(std::string info = "") : _info(std::move(info)) {}

  public:
    //! full constructor
    Exception(std::string info, std::string file, unsigned int line)
        : _info(std::move(info)), _file(std::move(file)), _line(line) {}

    /* ---------------------------------------------------------------------- */
    /*  Methods */
    /* ---------------------------------------------------------------------- */
  public:
    const char * what() const noexcept override { return _info.c_str(); }

    virtual std::string info() const noexcept {
      std::stringstream stream;
      stream << debug::demangle(typeid(*this).name()) << " : " << _info << " ["
             << _file << ":" << _line << "]";
      return stream.str();
    }

  public:
    void setInfo(const std::string & info) { _info = info; }
    void setFile(const std::string & file) { _file = file; }
    void setLine(unsigned int line) { _line = line; }
    void setModule(const std::string & module) { _module = module; }

    void setBacktrace(const std::vector<std::string> & backtrace) {
      backtrace_ = backtrace;
    }

    decltype(auto) backtrace() const { return backtrace_; }
    /* ---------------------------------------------------------------------- */
    /* Class Members                                                          */
    /* ---------------------------------------------------------------------- */
  protected:
    /// exception description and additionals
    std::string _info;

  private:
    /// file it is thrown from
    std::string _file;

    /// line it is thrown from
    unsigned int _line{0};

    /// module in which exception was raised
    std::string _module{"core"};

    std::vector<std::string> backtrace_;
  };

  class CriticalError : public Exception {};
  class AssertException : public Exception {};
  class NotImplementedException : public Exception {};

  /// standard output stream operator
  inline std::ostream & operator<<(std::ostream & stream,
                                   const Exception & _this) {
    stream << _this.what();
    return stream;
  }

  /* --------------------------------------------------------------------------
   */
  class Debugger {
  public:
    Debugger() noexcept;
    virtual ~Debugger();
    Debugger(const Debugger &) = default;
    Debugger & operator=(const Debugger &) = default;
    Debugger(Debugger &&) noexcept = default;
    Debugger & operator=(Debugger &&) noexcept = default;

    static void exit(int status) __attribute__((noreturn));

    void throwException(const std::string & info, const std::string & file,
                        unsigned int line, bool /*silent*/,
                        const std::string & /*location*/,
                        const std::string & module) const noexcept(false)
        __attribute__((noreturn));

    /*----------------------------------------------------------------------- */
    template <class Except>
    void throwCustomException(Except ex, const std::string & info,
                              const std::string & file, unsigned int line,
                              const std::string & module) const noexcept(false)
        __attribute__((noreturn));
    /*----------------------------------------------------------------------- */
    template <class Except>
    void throwCustomException(Except ex, const std::string & file,
                              unsigned int line,
                              const std::string & module_) const noexcept(false)
        __attribute__((noreturn));

    void printMessage(const std::string & prefix, const DebugLevel & level,
                      const std::string & info,
                      const std::string & module_) const;

    void setOutStream(std::ostream & out) { cout = &out; }
    std::ostream & getOutStream() { return *cout; }

  public:
    void setParallelContext(int rank, int size);

    void setDebugLevel(const DebugLevel & level);
    const DebugLevel & getDebugLevel() const;

    void setLogFile(const std::string & filename);
    std::ostream & getOutputStream();

    inline bool testLevel(const DebugLevel & level,
                          const std::string & module = "core") const {
      auto level_reached = (this->level >= (level));
      auto correct_module =
          (level <= dblCritical) or (modules_to_debug.empty()) or
          (modules_to_debug.find(module) != modules_to_debug.end());
      return level_reached and correct_module;
    }

    void printBacktrace(bool on_off) { this->print_backtrace = on_off; }
    bool printBacktrace() const { return this->print_backtrace; }

    void addModuleToDebug(const std::string & id) {
      modules_to_debug.insert(id);
    }
    void removeModuleToDebug(const std::string & id) {
      auto it = modules_to_debug.find(id);
      if (it != modules_to_debug.end()) {
        modules_to_debug.erase(it);
      }
    }

    void listModules() {
      for (const auto & module_ : modules_to_debug) {
        (*cout) << module_ << std::endl;
      }
    }

  private:
    std::string parallel_context;
    std::ostream * cout;
    bool file_open;
    DebugLevel level;
    bool print_backtrace;
    std::set<std::string> modules_to_debug;
  };

  extern Debugger debugger; // NOLINT
} // namespace debug

/* -------------------------------------------------------------------------- */
#define AKANTU_STRINGIZE_(str) #str
#define AKANTU_STRINGIZE(str) AKANTU_STRINGIZE_(str)
/* -------------------------------------------------------------------------- */
#define AKANTU_DEBUG_MODULE AKANTU_STRINGIZE(AKANTU_MODULE)
/* -------------------------------------------------------------------------- */
#define AKANTU_STRINGSTREAM_IN(_str, _sstr)                                    \
  ;                                                                            \
  do {                                                                         \
    std::stringstream _dbg_s_info;                                             \
    _dbg_s_info << _sstr; /* NOLINT */                                         \
    (_str) = _dbg_s_info.str();                                                \
  } while (false)

/* -------------------------------------------------------------------------- */
#define AKANTU_EXCEPTION(info) AKANTU_EXCEPTION_(info, false)

#define AKANTU_SILENT_EXCEPTION(info) AKANTU_EXCEPTION_(info, true)

#define AKANTU_EXCEPTION_(info, silent)                                        \
  do {                                                                         \
    std::stringstream _dbg_str;                                                \
    _dbg_str << info; /* NOLINT */                                             \
    std::stringstream _dbg_loc;                                                \
    _dbg_loc << AKANTU_LOCATION;                                               \
    ::akantu::debug::debugger.throwException(_dbg_str.str(), __FILE__,         \
                                             __LINE__, silent, _dbg_loc.str(), \
                                             AKANTU_DEBUG_MODULE);             \
  } while (false)

#define AKANTU_CUSTOM_EXCEPTION_INFO(ex, info)                                 \
  do {                                                                         \
    std::stringstream _dbg_str;                                                \
    _dbg_str << info; /* NOLINT */                                             \
    ::akantu::debug::debugger.throwCustomException(                            \
        ex, _dbg_str.str(), __FILE__, __LINE__, AKANTU_DEBUG_MODULE);          \
  } while (false)

#define AKANTU_CUSTOM_EXCEPTION(ex)                                            \
  do {                                                                         \
    ::akantu::debug::debugger.throwCustomException(ex, __FILE__, __LINE__,     \
                                                   AKANTU_DEBUG_MODULE);       \
  } while (false)

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_NDEBUG
#define AKANTU_DEBUG_TEST(level) (false)
#define AKANTU_DEBUG_LEVEL_IS_TEST()                                           \
  (::akantu::debug::debugger.testLevel(dblTest, AKANTU_DEBUG_MODULE))
#define AKANTU_DEBUG(level, info)
#define AKANTU_DEBUG_(pref, level, info)
#define AKANTU_DEBUG_IN()
#define AKANTU_DEBUG_OUT()
#define AKANTU_DEBUG_INFO(info)
#define AKANTU_DEBUG_WARNING(info)
#define AKANTU_DEBUG_TRACE(info)
#define AKANTU_DEBUG_ASSERT(test, info)
#define AKANTU_ERROR(info)                                                     \
  AKANTU_CUSTOM_EXCEPTION_INFO(::akantu::debug::CriticalError(), info)
/* -------------------------------------------------------------------------- */
#else
#define AKANTU_DEBUG(level, info) AKANTU_DEBUG_("   ", level, info)

#define AKANTU_DEBUG_(pref, level, info)                                       \
  do {                                                                         \
    std::string _dbg_str;                                                      \
    AKANTU_STRINGSTREAM_IN(_dbg_str,                                           \
                           info << " " << AKANTU_LOCATION); /* NOLINT */       \
    ::akantu::debug::debugger.printMessage(pref, level, _dbg_str,              \
                                           AKANTU_DEBUG_MODULE);               \
  } while (false)

#define AKANTU_DEBUG_TEST(level)                                               \
  (::akantu::debug::debugger.testLevel(level, AKANTU_DEBUG_MODULE))

#define AKANTU_DEBUG_LEVEL_IS_TEST()                                           \
  (::akantu::debug::debugger.testLevel(dblTest))

#define AKANTU_DEBUG_IN()                                                      \
  AKANTU_DEBUG_(                                                               \
      "==>", ::akantu::dblIn,                                                  \
      __func__                                                                 \
          << "()") // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,
                   // bugprone-lambda-function-name)

#define AKANTU_DEBUG_OUT()                                                     \
  AKANTU_DEBUG_(                                                               \
      "<==", ::akantu::dblOut,                                                 \
      __func__                                                                 \
          << "()") // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,
                   // bugprone-lambda-function-name)

#define AKANTU_DEBUG_INFO(info) AKANTU_DEBUG_("---", ::akantu::dblInfo, info)

#define AKANTU_DEBUG_WARNING(info)                                             \
  AKANTU_DEBUG_("/!\\", ::akantu::dblWarning, info)

#define AKANTU_DEBUG_TRACE(info) AKANTU_DEBUG_(">>>", ::akantu::dblTrace, info)

#define AKANTU_DEBUG_ASSERT(test, info)                                        \
  do {                                                                         \
    if (not(test))                                                             \
      AKANTU_CUSTOM_EXCEPTION_INFO(::akantu::debug::AssertException(),         \
                                   "assert [" << #test << "] "                 \
                                              << info); /* NOLINT */           \
  } while (false)

#define AKANTU_ERROR(info)                                                     \
  do {                                                                         \
    AKANTU_DEBUG_("!!! ", ::akantu::dblError, info);                           \
    AKANTU_CUSTOM_EXCEPTION_INFO(::akantu::debug::CriticalError(),             \
                                 info); /* NOLINT */                           \
  } while (false)
#endif // AKANTU_NDEBUG

#define AKANTU_TO_IMPLEMENT()                                                  \
  AKANTU_CUSTOM_EXCEPTION_INFO(                                                \
      ::akantu::debug::NotImplementedException(),                              \
      __func__                                                                 \
          << " : not implemented yet !") // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,
                                         // bugprone-lambda-function-name)

/* -------------------------------------------------------------------------- */

namespace debug {
  /* ------------------------------------------------------------------------ */
  template <class Except>
  void
  Debugger::throwCustomException(Except ex, const std::string & info,
                                 const std::string & file, unsigned int line,
                                 const std::string & module_) const
      noexcept(false) {
    ex.setInfo(info);
    ex.setFile(file);
    ex.setLine(line);
    ex.setModule(module_);
    if (::akantu::debug::debugger.printBacktrace()) {
      ex.setBacktrace(::akantu::debug::getBacktrace());
    }
    throw ex;
  }
  /* ------------------------------------------------------------------------ */
  template <class Except>
  void Debugger::throwCustomException(Except ex, const std::string & file,
                                      unsigned int line,
                                      const std::string & module_) const
      noexcept(false) {
    ex.setFile(file);
    ex.setLine(line);
    ex.setModule(module_);
    if (::akantu::debug::debugger.printBacktrace()) {
      ex.setBacktrace(::akantu::debug::getBacktrace());
    }
    throw ex;
  }
} // namespace debug
} // namespace akantu

#endif /* AKANTU_ERROR_HH_ */
