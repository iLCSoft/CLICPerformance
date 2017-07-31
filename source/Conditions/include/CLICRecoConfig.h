#ifndef CLICRecoConfig_h
#define CLICRecoConfig_h 1

#include "marlin/Processor.h"
#include "marlin/EventModifier.h"

class CLICRecoConfig : public marlin::Processor, public marlin::EventModifier {

public:

  virtual marlin::Processor*  newProcessor() { return new CLICRecoConfig; }
  virtual std::string const& name() const { return marlin::Processor::name(); };

  CLICRecoConfig() ;
  CLICRecoConfig(const CLICRecoConfig&) = delete;
  CLICRecoConfig& operator=(const CLICRecoConfig&) = delete;

  // Initialisation - run at the beginning to start histograms, etc.
  virtual void init() ;

  // Called at the beginning of every run
  virtual void processRunHeader( LCRunHeader* run ) ;

  // Run over each event - the main algorithm
  virtual void modifyEvent( LCEvent* evt ) ;
  virtual void processEvent( LCEvent* evt ) ;

  // Run at the end of each event
  virtual void check( LCEvent* evt ) ;

  // Called at the very end for cleanup, histogram saving, etc.
  virtual void end() ;

protected:

  // Parameters
  const std::vector< std::string > m_trackingPossibleOptions{ "Truth", "Conformal" };
  std::string m_trackingOption="Truth";

  // Parameters
  std::vector< std::string > m_overlayPossibleOptions{ "False", "350GeV", "380GeV", "420GeV", "500GeV", "1.4TeV", "3TeV" };
  std::map<std::string, bool> m_overlay{};
  std::string m_overlayChoice="False";

  bool m_truthTracking = true;
  bool m_conformalTracking = false;


protected:
  void checkOptions( std::string const& optionValue, std::vector<std::string> const& options ) const;

  template< class T > static void toString( std::stringstream& stream, std::vector<T> const& vector) {
    unsigned int counter = 0;
    for (auto const& entry: vector) {
      stream << entry;
      if( counter < vector.size() - 1 ) {
        stream << ", ";
      }
      counter++;
    }
  }

} ;

#endif
