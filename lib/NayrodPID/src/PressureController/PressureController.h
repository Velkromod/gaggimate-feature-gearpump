// PressureController.h
#ifndef PRESSURE_CONTROLLER_H
#define PRESSURE_CONTROLLER_H

#ifndef M_PI
static constexpr float M_PI = 3.14159265358979323846f;
#endif

#include "SimpleKalmanFilter/SimpleKalmanFilter.h"
#include <algorithm>
#include <cmath>

class PressureController {
  private:
    static void applyLowPassFilter(float *filteredValue, float rawValue, float cutoffFreq, float dt);

  public:
    enum class ControlMode { POWER, PRESSURE, FLOW };

    PressureController(float dt, float *_rawPressureSetpoint, float *_rawFlowSetpoint, float *sensorOutput,
                       float *controllerOutput, int *valveStatus);

    void initSetpointFilter(float val = 0.0f);

    void setFlowLimit(float lim) { /* Flow limit not currently implemented */ }
    void setPressureLimit(float lim) { /* Pressure limit not currently implemented */ }

    void update(ControlMode mode);
    void tare();
    void reset();

    float getCoffeeOutputEstimate() { return std::fmax(0.0f, _coffeeOutput); }
    void setPumpFlowCoeff(float oneBarFlow, float nineBarFlow);
    void setPumpFlowPolyCoeffs(float a, float b, float c, float d);
    float getPumpFlowRate() { return exportPumpFlowRate; }
    float getCoffeeFlowRate() { return *_valveStatus == 1 ? _coffeeFlowRate : 0.0f; }
    float getPuckResistance() { return _puckResistance; }

    float getRawPressure() const { return _rawPressure ? *_rawPressure : 0.0f; }
    float getFilteredPressure() const { return _filteredPressureSensor; }

    float getRawPressureSetpoint() const { return _rawPressureSetpoint ? *_rawPressureSetpoint : 0.0f; }
    float getFilteredPressureSetpoint() const { return _filteredSetpoint; }
    float getFilteredSetpointDerivative() const { return _filteredSetpointDerivative; }

    float getControlOutput() const { return _ctrlOutput ? *_ctrlOutput : 0.0f; }
    float getFilteredPressureDerivative() const { return _filteredPressureDerivative; }
    float getPumpDutyCycleInternal() const { return _pumpDutyCycle; }

    float getRawPressureControlOutput() const { return _lastRawPressureOutput; }
    float getAppliedPressureControlOutput() const { return _lastAppliedPressureOutput; }
    float getPressureLimiterDelta() const { return _lastRawPressureOutput - _lastAppliedPressureOutput; }
    bool isPressureLimiterActiveUp() const { return _pressureLimiterActiveUp; }
    bool isPressureLimiterActiveDown() const { return _pressureLimiterActiveDown; }

    float getUnclampedPressureControlOutput() const { return _lastUnclampedPressureOutput; }
    float getClampedPressureControlOutput() const { return _lastClampedPressureOutput; }
    float getPressureClampDelta() const { return _lastUnclampedPressureOutput - _lastClampedPressureOutput; }
    bool isPressureClampActiveHigh() const { return _pressureClampActiveHigh; }
    bool isPressureClampActiveLow() const { return _pressureClampActiveLow; }

    float getPressureFeedbackOutput() const { return _lastFeedbackOutput; }
    float getPressureFeedforwardHoldOutput() const { return _lastFeedforwardHoldOutput; }
    float getPressureFeedforwardDynamicOutput() const { return _lastFeedforwardDynamicOutput; }
    float getPressureFeedforwardTotalOutput() const { return _lastFeedforwardTotalOutput; }

    float getErrorIntegral() const { return _errorIntegral; }

    void setDeadVolume(float deadVol) { _puckSaturatedVolume = deadVol; }

  private:
    float getPumpDutyCycleForPressure();
    void virtualScale();
    void filterSensor();
    void filterSetpoint(float rawSetpoint);

    float pumpFlowModel(float alpha = 100.0f) const;
    float getAvailableFlow() const;
    float getAvailableFlowAtPressure(float pressure) const;
    float getPumpDutyCycleForFlowRate() const;

    float getScheduledSetpointWeight(float pressureRef) const;
    float getPressureFeedforwardHold(float pressureRef) const;
    float getPressureFeedforwardDynamicUnfiltered(float pressureRef, float pressureRefDerivative, float absError) const;

    float _dt = 1.0f;

    float *_rawPressureSetpoint = nullptr;
    float *_rawFlowSetpoint = nullptr;
    float *_rawPressure = nullptr;
    float *_ctrlOutput = nullptr;
    int *_valveStatus = nullptr;

    float _filteredPressureSensor = 0.0f;
    float _filteredSetpoint = 0.0f;
    float _filteredSetpointDerivative = 0.0f;
    float _filteredPressureDerivative = 0.0f;

    float _setpointFilterFreq = 1.6f;
    float _setpointFilterDamping = 1.2f;
    bool _setpointFilterInitialized = false;

    const float _systemCompliance = 1.4f;
    float _puckResistance = 1e7f;
    const float _maxPressure = 15.0f;
    const float _maxPressureRate = 13.0f;
    float _pumpFlowCoefficients[4] = {0.0f, 0.0f, -0.5854f, 10.79f};

    float _commutationGain = 0.180f;
    float _convergenceGain = 1.03f;
    float _epsilonCoefficient = 0.40f;
    float _deadbandCoefficient = 0.06f;
    float _integralGain = 0.08f;

    float _setpointWeightLowPressure = 1.0f;
    float _setpointWeightHighPressure = 1.0f;
    float _setpointWeightScheduleStartBar = 3.0f;
    float _setpointWeightScheduleEndBar = 9.0f;

    float _feedforwardBiasPct = 0.0f;
    float _feedforwardPressureGainPctPerBar = 0.0f;
    float _feedforwardHoldMaxPct = 0.0f;

    float _feedforwardRampGain = 0.28f;
    float _feedforwardDynamicMaxPct = 6.0f;
    float _feedforwardDynamicFilterFreq = 1.0f;
    float _filteredFeedforwardDynamicOutput = 0.0f;
    float _lastFeedforwardDynamicAppliedOutput = 0.0f;

    float _feedforwardDynamicRiseRate = 50.0f;
    float _feedforwardDynamicDropRate = 70.0f;

    float _ffPressureGateStartBar = 1.5f;
    float _ffPressureGateFullBar = 4.0f;
    float _ffPressureTaperStartBar = 7.6f;
    float _ffPressureTaperEndBar = 9.2f;
    float _ffErrorGateStartBar = 0.25f;
    float _ffErrorGateFullBar = 0.90f;

    float _previousPressure = 0.0f;
    float _errorIntegral = 0.0f;
    float _pumpDutyCycle = 0.0f;

    float _lastFeedbackOutput = 0.0f;
    float _lastFeedforwardHoldOutput = 0.0f;
    float _lastFeedforwardDynamicOutput = 0.0f;
    float _lastFeedforwardTotalOutput = 0.0f;

    float _lastUnclampedPressureOutput = 0.0f;
    float _lastClampedPressureOutput = 0.0f;

    float _lastRawPressureOutput = 0.0f;
    float _lastAppliedPressureOutput = 0.0f;

    bool _pressureClampActiveHigh = false;
    bool _pressureClampActiveLow = false;
    bool _pressureLimiterActiveUp = false;
    bool _pressureLimiterActiveDown = false;

    const float _maxPressureOutputRiseRate = 190.0f;
    const float _maxPressureOutputDropRate = 255.0f;

    float _waterThroughPuckFlowRate = 0.0f;
    float _pumpFlowRate = 0.0f;
    float _pumpVolume = 0.0f;
    float _coffeeOutput = 0.0f;
    float _coffeeFlowRate = 0.0f;
    float _lastFilteredPressure = 0.0f;
    float _filterEstimatorFrequency = 1.0f;
    float _pressureFilterEstimator = 0.0f;
    float _puckSaturationVolume = 0.0f;
    float _puckSaturatedVolume = 45.0f;
    float _lastPuckConductance = 0.0f;
    float _puckConductance = 0.0f;
    float _puckConductanceDerivative = 0.0f;
    bool _puckState[3] = {};
    int _puckCounter = 0;

    float exportPumpFlowRate = 0.0f;
    SimpleKalmanFilter *_pressureKalmanFilter = nullptr;
};

#endif // PRESSURE_CONTROLLER_H
