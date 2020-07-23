#ifndef EVENT_SHAPE_H
#define EVENT_SHAPE_H

// STL headers
#include <algorithm>
#include <cmath>

#include "SampleAnalyzer/Commons/DataFormat/EventFormat.h"
#include "SampleAnalyzer/Commons/DataFormat/SampleFormat.h"
#include "SampleAnalyzer/Commons/Vector/MAVector3.h"

namespace MA5
{
  class EventShape
  {
    private:

        MAdouble64 Sphericity_;
        MAdouble64 Aplanarity_;
        MAdouble64 TransverseSphericity_;
        
        MAdouble64 thrust_;
        MAdouble64 thrustMajor_;
        MAdouble64 thrustMinor_;

        std::vector<MAVector3> thrustAxes_;

    public :
        /// Constructor without argument
        EventShape()
        {
            Sphericity_           = -1;
            Aplanarity_           = -1;
            TransverseSphericity_ = -1;
            thrust_               = -1;
            thrustMajor_          = -1;
            thrustMinor_          = -1;
            thrustAxes_.clear();
        }
        
        /// Destructor
        ~EventShape() {}
        
        void Reset()
        {
            Sphericity_           = -1;
            Aplanarity_           = -1;
            TransverseSphericity_ = -1;
            thrust_               = -1;
            thrustMajor_          = -1;
            thrustMinor_          = -1;
            thrustAxes_.clear();
        }

        void calculateSphericity(std::vector<const RecJetFormat*> Jets);
        void calculateThrust(std::vector<const RecJetFormat*> Jets);

        MAdouble64 Sphericity()  {return Sphericity_;}
        MAdouble64 Aplanarity()  {return Aplanarity_;}
        MAdouble64 TSphericity() {return TransverseSphericity_;}

        // The thrust axis
        MAdouble64 thrust() { return thrust_;}
        /// The thrust major scalar, \f$ M \f$, (thrust along thrust major axis).
        MAdouble64 thrustMajor() { return thrustMajor_; }
        /// The thrust minor scalar, \f$ m \f$, (thrust along thrust minor axis).
        MAdouble64 thrustMinor() { return thrustMinor_; }
        /// The oblateness, \f$ O = M - m \f$ .
        MAdouble64 oblateness() { return thrustMajor_ - thrustMinor_; }

        /// The thrust axis.
        MAVector3 thrustAxis() { return thrustAxes_[0]; }
        /// The thrust major axis (axis of max thrust perpendicular to thrust axis).
        MAVector3 thrustMajorAxis() { return thrustAxes_[1]; }
        /// The thrust minor axis (axis perpendicular to thrust and thrust major).
        MAVector3 thrustMinorAxis() { return thrustAxes_[2]; }

  };
}
#endif