#ifndef CLICRecoConfig_h
#define CLICRecoConfig_h 1

#include "marlin/Processor.h"
#include "marlin/EventModifier.h"

class CLICRecoConfig : public marlin::Processor, public marlin::EventModifier {

  using Choices = std::vector<std::string>;

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

  std::map<std::string, bool> m_options{};

  // Tracking
  const Choices m_trackingPossibleOptions{"Truth", "ConformalPlusExtrapolator", "Conformal"};
  std::string m_trackingChoice="Truth";

  // Overlay
  const Choices m_overlayPossibleOptions{"False", "91GeV", "350GeV", "365GeV",  "380GeV", "420GeV", "500GeV", "1.4TeV", "3TeV"};
  std::string m_overlayChoice="False";

  // BeamCal
  const Choices m_beamcalPossibleOptions{"3TeV", "380GeV"};
  std::string m_beamcalChoice="3TeV";

  std::vector<std::tuple<std::string, std::string, Choices>> m_allOptions =
    {{"Tracking", m_trackingChoice, m_trackingPossibleOptions},
     {"BeamCal",  m_beamcalChoice, m_beamcalPossibleOptions},
     {"Overlay",  m_overlayChoice, m_overlayPossibleOptions}};


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
