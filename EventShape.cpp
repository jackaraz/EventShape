#include "SampleAnalyzer/User/Analyzer/EventShape.h"

using namespace MA5;

void EventShape::calculateSphericity(std::vector<const RecJetFormat*> Jets)
{
    // Construction of the sphericity tensor, calculation of aplanarity
    // using Cardano algorithm
    if (Jets.size() == 0) return;
    MAdouble64 S12 = 0., S31 = 0., S23 = 0., S11 = 0., S22 = 0., S33 = 0., Stot = 0.;
    for (MAuint32 i = 0; i < Jets.size(); i++)
    {
        S11  += Jets[i]->px() * Jets[i]->px();
        S12  += Jets[i]->px() * Jets[i]->py();
        S22  += Jets[i]->py() * Jets[i]->py();
        S23  += Jets[i]->py() * Jets[i]->pz();
        S31  += Jets[i]->pz() * Jets[i]->px();
        S33  += Jets[i]->pz() * Jets[i]->pz();
        Stot += Jets[i]->pt() * Jets[i]->pt();
    }
    S11=S11/Stot; S12=S12/Stot; S22=S22/Stot; S23=S23/Stot; S31=S31/Stot; S33=S33/Stot;
    MAdouble64 Sii    = S11+S22+S33;
    MAdouble64 C0     = S11*S23*S23 + S22*S31*S31 + S33*S12*S12 - S11*S22*S33 - 2.*S31*S12*S23;
    MAdouble64 C1     = S11*S22 + S22*S33 + S11*S33 - S12*S12 - S23*S23 - S31*S31;
    MAdouble64 P      = Sii*Sii - 3.*C1;
    MAdouble64 Q      = Sii*(P-1.5*C1) - 13.5*C0;
    MAdouble64 phi    = atan2(sqrt(fabs(27.*(C1*C1/4.*(P-C1) + C0*(Q+6.75*C0)))),Q)/3.;
    MAdouble64 cth    = sqrt(std::fabs(P))*cos(phi);
    MAdouble64 sth    = sqrt(std::fabs(P))*sin(phi)/sqrt(3.);
    MAdouble64 a1     = (Sii-cth)/3.+sth;
    MAdouble64 a2     = (Sii-cth)/3.-sth;
    MAdouble64 a3     = (Sii-cth)/3.+cth;
    std::vector<MAdouble64> eig;
    eig.push_back(a1);
    eig.push_back(a2);
    eig.push_back(a3);
    std::sort(eig.begin(),eig.end());
    MAdouble64 lam3   = eig[0];
    MAdouble64 lam2   = eig[1];
    MAdouble64 lam1   = eig[2];

    MAdouble64 lamtot = lam1 + lam2 + lam3;
    lam1 = lam1/lamtot; lam2 = lam2/lamtot; lam3 = lam3/lamtot;

    Aplanarity_           = 1.5*lam3;
    Sphericity_           = 1.5*(lam2+lam3);
    TransverseSphericity_ = 2.0*lam2/(lam1+lam2);
}


MAbool sort_by_mag2(MAVector3 & a, MAVector3 & b) {return a.Mag2() < b.Mag2();}

void calcT(std::vector<MAVector3>& momenta, MAdouble64& val, MAVector3& axis)
{
    // This function implements the iterative algorithm as described in the
    // Pythia manual. We take eight (four) different starting vectors
    // constructed from the four (three) leading particles to make sure that
    // we don't find a local maximum.
    std::vector<MAVector3> p = momenta;
    std::sort(p.begin(), p.end(), sort_by_mag2);
    MAuint32 n = 3;
    std::vector< MAVector3 >  tmp_vec;
    std::vector< MAdouble64 > tmp_val;
    for (int i = 0 ; i < int(pow(2, n-1)); i++)
    {
        // Create an initial vector from the leading four jets
        MAVector3 foo;
        int sign = i;
        for (MAuint32 k = 0 ; k < n ; ++k)
        {
            (sign % 2) == 1 ? foo += p[k] : foo -= p[k];
            sign /= 2;
        }
        foo=foo.Unit();

        // Iterate
        MAdouble64 diff=999.;
        while (diff>1e-5)
        {
            MAVector3 foobar;
            for (MAuint32 k=0 ; k<p.size() ; k++)
                foo.Dot(p[k])>0 ? foobar+=p[k] : foobar-=p[k];
            diff=(foo-foobar.Unit()).Mag();
            foo=foobar.Unit();
        }

        MAdouble64 tmp = 0.;
        // Calculate the thrust value for the vector we found
        for (MAuint32 k=0 ; k<p.size() ; k++)
            tmp += std::fabs(foo.Dot(p[k]));

        // Store everything
        tmp_val.push_back(tmp);
        tmp_vec.push_back(foo);
    }

    // Pick the solution with the largest thrust
    for (MAuint32 i=0 ; i<tmp_vec.size() ; i++)
    {
        if (tmp_val[i] > val)
        {
            val  = tmp_val[i];
            axis = tmp_vec[i];
        }
    }
}




void EventShape::calculateThrust(std::vector<const RecJetFormat*> Jets)
{
    if (Jets.size()<2)
    {
        for (MAuint32 i=0; i<3; i++)
            thrustAxes_.push_back(MAVector3(0.,0.,0.));
    }
    else if (Jets.size() == 2)
    {
        MAVector3 tmp_axis = MAVector3(0.,0.,0.);
        thrust_      = 1.0;
        thrustMajor_ = 0.0;
        thrustMinor_ = 0.0;
        tmp_axis = Jets[0]->momentum().Vect().Unit();
        if (tmp_axis.Z()<0) tmp_axis = -tmp_axis;
        thrustAxes_.push_back(tmp_axis);
        if (tmp_axis.Z()<0.75)
        {    thrustAxes_.push_back((tmp_axis.Cross(MAVector3(0.,0.,1.))).Unit());}
        else
        {    thrustAxes_.push_back((tmp_axis.Cross(MAVector3(0.,1.,0.))).Unit());}
        thrustAxes_.push_back(thrustAxes_[0].Cross(thrustAxes_[1]));
    }
    else
    {
        std::vector<MAVector3> momenta;
        MAdouble64 MomentumSum;
        for (MAuint32 i=0; i<Jets.size(); i++)
        {
            momenta.push_back(Jets[i]->momentum().Vect());
            MomentumSum += Jets[i]->momentum().Vect().Mag();
        }

        // Get Thrust
        MAdouble64 val = 0.;
        MAVector3  axis;
        calcT(momenta, val, axis);
        thrust_ = val / MomentumSum;
        if (axis.Z() < 0) axis = -axis;
        axis = axis.Unit();
        thrustAxes_.push_back(axis);

        // Get thrust major
        std::vector<MAVector3> threeMomenta;
        // Get the part of each 3-momentum which is perpendicular 
        // to the thrust axis
        for (MAuint32 i=0; i<Jets.size(); i++)
        {
            MAVector3 current_jet_mom = Jets[i]->momentum().Vect();
            MAVector3 vpar = current_jet_mom.Dot(axis.Unit()) * axis.Unit();
            threeMomenta.push_back(current_jet_mom - vpar);
        }
        val = 0.;
        axis = MAVector3(0.,0.,0.);
        calcT(threeMomenta, val, axis);
        thrustMajor_ = val / MomentumSum;
        if (axis.X() < 0) axis = -axis;
        axis = axis.Unit();
        thrustAxes_.push_back(axis);

        // Get thrust minor
        if (thrustAxes_[0].Dot(thrustAxes_[1]) < 1e-10)
        {
            axis = thrustAxes_[0].Cross(thrustAxes_[1]);
            thrustAxes_.push_back(axis);
            val = 0.0;
            for (MAuint32 i=0; i<Jets.size(); i++)
            {
                MAVector3 current_jet_mom = Jets[i]->momentum().Vect();
                val += std::fabs(current_jet_mom.Dot(axis));
            }
            thrustMinor_ = val / MomentumSum;
        }
        else
        {
            thrustMinor_ = -1.0;
            thrustAxes_.push_back(MAVector3(0,0,0));
        }
    }
}














