#include "PressureController.h"
#include "SimpleKalmanFilter/SimpleKalmanFilter.h"
#include "esp_log.h"

#include <algorithm>
#include <math.h>

// Generic first-order low-pass helper reused by estimator / telemetry / FF path.
void PressureController::applyLowPassFilter(float *filteredValue, float rawValue, float cutoffFreq, float dt) {
    if (filteredValue == nullptr)
        return;

    float alpha = dt / (1.0f / (2.0f * M_PI * cutoffFreq) + dt);
    *filteredValue = alpha * rawValue + (1.0f - alpha) * (*filteredValue);
}

PressureController::PressureController(float dt, float *rawPressureSetpoint, float *rawFlowSetpoint, float *sensorOutput,
                                       float *controllerOutput, int *valveStatus) {
    this->_rawPressureSetpoint = rawPressureSetpoint;
    this->_rawFlowSetpoint = rawFlowSetpoint;
    this->_rawPressure = sensorOutput;
    this->_ctrlOutput = controllerOutput;
    this->_valveStatus = valveStatus;
    this->_dt = dt;

    // Kalman process noise term previously retuned for this gear-pump branch.
    this->_pressureKalmanFilter = new SimpleKalmanFilter(0.1f, 10.0f, powf(8 * _dt, 2));
    this->_previousPressure = *sensorOutput;
}

void PressureController::filterSetpoint(float rawSetpoint) {
    if (!_setpointFilterInitialized)
        initSetpointFilter();

    float omega = 2.0f * M_PI * _setpointFilterFreq;
    float d2r =
        (omega * omega) * (rawSetpoint - _filteredSetpoint) - 2.0f * _setpointFilterDamping * omega * _filteredSetpointDerivative;

    _filteredSetpointDerivative += std::clamp(d2r * _dt, -_maxPressureRate, _maxPressureRate);
    _filteredSetpoint += _filteredSetpointDerivative * _dt;
}

void PressureController::initSetpointFilter(float val) {
    _filteredSetpoint = *_rawPressureSetpoint;
    if (val != 0.0f)
        _filteredSetpoint = val;

    _filteredSetpointDerivative = 0.0f;
    _setpointFilterInitialized = true;
}

void PressureController::filterSensor() {
    float newFiltered = this->_pressureKalmanFilter->updateEstimate(*_rawPressure);

    float pressureDerivative = (newFiltered - _lastFilteredPressure) / _dt;
    applyLowPassFilter(&_filteredPressureDerivative, pressureDerivative, _filterEstimatorFrequency, _dt);

    _lastFilteredPressure = newFiltered;
    _filteredPressureSensor = newFiltered;
}

void PressureController::update(ControlMode mode) {
    filterSetpoint(*_rawPressureSetpoint);
    filterSensor();

    if ((mode == ControlMode::FLOW || mode == ControlMode::PRESSURE) && *_rawPressureSetpoint > 0.0f &&
        *_rawFlowSetpoint > 0.0f) {
        float flowOutput = getPumpDutyCycleForFlowRate();
        float pressureOutput = getPumpDutyCycleForPressure();
        *_ctrlOutput = std::min(flowOutput, pressureOutput);

        if (flowOutput < pressureOutput) {
            // When the flow branch dominates, clear the pressure integrator so
            // the pressure loop does not keep hidden accumulated error.
            _errorIntegral = 0.0f;
        }
    } else if (mode == ControlMode::FLOW) {
        *_ctrlOutput = getPumpDutyCycleForFlowRate();
    } else if (mode == ControlMode::PRESSURE) {
        *_ctrlOutput = getPumpDutyCycleForPressure();
    }

    virtualScale();
}

float PressureController::pumpFlowModel(float alpha) const {
    const float availableFlow = getAvailableFlow();
    return availableFlow * alpha / 100.0f;
}

// Helper that evaluates the pump flow model at an arbitrary pressure.
// This is used by the dynamic FF path to convert desired dP/dt into an
// estimated actuator contribution based on the local operating point.
float PressureController::getAvailableFlowAtPressure(float pressure) const {
    const float P = std::max(0.0f, pressure);
    const float P2 = P * P;
    const float P3 = P2 * P;

    const float Q =
        _pumpFlowCoefficients[0] * P3 + _pumpFlowCoefficients[1] * P2 + _pumpFlowCoefficients[2] * P + _pumpFlowCoefficients[3];

    return Q;
}

float PressureController::getAvailableFlow() const {
    return getAvailableFlowAtPressure(_filteredPressureSensor);
}

float PressureController::getPumpDutyCycleForFlowRate() const {
    const float availableFlow = getAvailableFlow();
    if (availableFlow <= 0.0f) {
        return 0.0f;
    }

    float duty = (*_rawFlowSetpoint / availableFlow) * 100.0f;
    return std::clamp(duty, 0.0f, 100.0f);
}

void PressureController::setPumpFlowCoeff(float oneBarFlow, float nineBarFlow) {
    _pumpFlowCoefficients[0] = 0.0f;
    _pumpFlowCoefficients[1] = 0.0f;
    _pumpFlowCoefficients[2] = (nineBarFlow - oneBarFlow) / 8.0f;
    _pumpFlowCoefficients[3] = oneBarFlow - _pumpFlowCoefficients[2] * 1.0f;
}

void PressureController::setPumpFlowPolyCoeffs(float a, float b, float c, float d) {
    _pumpFlowCoefficients[0] = a;
    _pumpFlowCoefficients[1] = b;
    _pumpFlowCoefficients[2] = c;
    _pumpFlowCoefficients[3] = d;
}

void PressureController::tare() {
    _coffeeOutput = 0.0f;
    _pumpVolume = 0.0f;
    _puckSaturationVolume = 0.0f;
    _puckState[0] = false;
    _puckState[1] = false;
    _puckState[2] = false;
    _puckCounter = 0;
    _pumpFlowRate = 0.0f;
    _puckConductanceDerivative = 0.0f;
    _coffeeFlowRate = 0.0f;
    _puckResistance = INFINITY;
}

void PressureController::virtualScale() {
    float newPumpFlowRate = pumpFlowModel(*_ctrlOutput);
    applyLowPassFilter(&_pumpFlowRate, newPumpFlowRate, _filterEstimatorFrequency, _dt);
    _pumpVolume += _pumpFlowRate * _dt;
    applyLowPassFilter(&exportPumpFlowRate, newPumpFlowRate, _filterEstimatorFrequency / 2.0f, _dt);

    float effectiveCompliance = 3.0f / fmax(0.2f, _filteredPressureSensor); // ml*s/bar
    float flowRaw = _pumpFlowRate - effectiveCompliance * _filteredPressureDerivative;

    applyLowPassFilter(&_waterThroughPuckFlowRate, flowRaw, 0.3f, _dt);
    if (_waterThroughPuckFlowRate > 0.0f && *_valveStatus == 1 && _filteredPressureSensor > 0.8f) {
        _puckCounter++;
        _puckSaturationVolume += _waterThroughPuckFlowRate * _dt;

        _puckConductance = _waterThroughPuckFlowRate / sqrtf(_filteredPressureSensor);

        if (_puckCounter <= 1)
            _lastPuckConductance = _puckConductance;

        float newPuckConductanceDerivative = (_puckConductance - _lastPuckConductance) / _dt;
        applyLowPassFilter(&_puckConductanceDerivative, newPuckConductanceDerivative, 0.2f, _dt);
        _lastPuckConductance = _puckConductance;

        if (_puckConductanceDerivative < -0.5f && float(_puckCounter) * _dt > 1.0f)
            _puckState[0] = true;

        int timeStamp = 0;
        if (_puckState[0] && _puckConductanceDerivative > -0.1f && !_puckState[1]) {
            _puckState[1] = true;
            timeStamp = _puckCounter;
        }

        _puckResistance = 1.0f / _puckConductance;
        if (_puckState[1]) {
            if (_puckCounter == timeStamp) {
                _coffeeFlowRate = _waterThroughPuckFlowRate;
                _puckResistance = 1.0f / _puckConductance;
            }

            applyLowPassFilter(&_puckResistance, 1.0f / _puckConductance, 0.1f, _dt);
            applyLowPassFilter(&_coffeeFlowRate, _waterThroughPuckFlowRate, 0.2f, _dt);

            if (!_puckState[2]) {
                float timeMissedDrops = 2.0f;
                float missedDrops = _coffeeFlowRate * timeMissedDrops / 2.0f;
                _puckState[2] = true;
                _coffeeOutput += _coffeeFlowRate * _dt + missedDrops;
            } else {
                _coffeeOutput += _coffeeFlowRate * _dt;
            }
        }
    }
}

// 2-DOF weighting hook kept in place, but this iteration still keeps weighting disabled.
// The goal is to isolate the effect of dynamic FF changes.
float PressureController::getScheduledSetpointWeight(float pressureRef) const {
    (void)pressureRef;
    return 1.0f;
}

// Static hold FF remains disabled.
float PressureController::getPressureFeedforwardHold(float pressureRef) const {
    (void)pressureRef;
    return 0.0f;
}

// Continuous dynamic feedforward.
// This path:
// 1) acts only during positive setpoint ramps
// 2) uses the local pump-flow model and system compliance
// 3) applies a smooth pressure window
// 4) stays active while pressure is still below target
// 5) adds a bounded gamma boost during genuine ramps to prevent premature feedback unloading
float PressureController::getPressureFeedforwardDynamicUnfiltered(float pressureRef, float pressureRefDerivative, float signedError) const {
    // Only assist during upward setpoint motion.
    float risingRefDerivative = std::max(0.0f, pressureRefDerivative);
    if (risingRefDerivative <= 0.0f) {
        return 0.0f;
    }

    float availableFlowAtRef = fmaxf(getAvailableFlowAtPressure(std::max(0.0f, pressureRef)), 1e-3f);

    // Compliance-based required flow:
    // Q_required = C_eq * dP_ref/dt
    float requiredPercent = 100.0f * (_systemCompliance * risingRefDerivative) / availableFlowAtRef;

    auto smoothstep = [](float edge0, float edge1, float x) -> float {
        if (edge1 <= edge0) {
            return (x >= edge1) ? 1.0f : 0.0f;
        }
        float t = (x - edge0) / (edge1 - edge0);
        t = std::clamp(t, 0.0f, 1.0f);
        return t * t * (3.0f - 2.0f * t);
    };

    // Pressure window:
    // fade in through low pressure, stay active in the mid-ramp, fade out near the top.
    float pressureGateIn = smoothstep(_ffPressureGateStartBar, _ffPressureGateFullBar, pressureRef);
    float pressureGateOut = 1.0f - smoothstep(_ffPressureTaperStartBar, _ffPressureTaperEndBar, pressureRef);
    float pressureWeight = std::clamp(pressureGateIn * pressureGateOut, 0.0f, 1.0f);

    // Above-target shutoff:
    // signedError = P - Pref.
    // Keep FF active while pressure is below or near target, and fade it out only when pressure rises above target.
    float aboveWeight = 1.0f - smoothstep(_ffAboveGateStartBar, _ffAboveGateFullBar, signedError);

    // Ramp-dependent gamma boost:
    // compensate the natural unloading of the feedback branch during a real positive ramp.
    float rampWeight = smoothstep(_ffRampDerivGateStart, _ffRampDerivGateFull, risingRefDerivative);
    float gamma = 1.0f + _ffGammaBoostMax * rampWeight * pressureWeight;

    float ff = gamma * _feedforwardRampGain * requiredPercent * pressureWeight * aboveWeight;
    return std::clamp(ff, 0.0f, _feedforwardDynamicMaxPct);
}

float PressureController::getPumpDutyCycleForPressure() {
    if (*_rawPressureSetpoint < 0.2f) {
        initSetpointFilter();
        _errorIntegral = 0.0f;
        _pumpDutyCycle = 0.0f;

        if (_ctrlOutput != nullptr) {
            *_ctrlOutput = 0.0f;
        }

        // Reset all telemetry / FF / limiter state when pressure control is inactive.
        _lastFeedbackOutput = 0.0f;
        _lastFeedforwardHoldOutput = 0.0f;
        _lastFeedforwardDynamicOutput = 0.0f;
        _lastFeedforwardTotalOutput = 0.0f;
        _filteredFeedforwardDynamicOutput = 0.0f;
        _lastFeedforwardDynamicAppliedOutput = 0.0f;

        _lastUnclampedPressureOutput = 0.0f;
        _lastClampedPressureOutput = 0.0f;
        _lastRawPressureOutput = 0.0f;
        _lastAppliedPressureOutput = 0.0f;

        _pressureClampActiveHigh = false;
        _pressureClampActiveLow = false;
        _pressureLimiterActiveUp = false;
        _pressureLimiterActiveDown = false;

        _previousPressure = 0.0f;
        return 0.0f;
    }

    float P = _filteredPressureSensor;
    float P_ref = _filteredSetpoint;

    // Full error is preserved for the integral term.
    float error = P - P_ref;

    // Weighted path kept available, but still neutralized in this iteration.
    float beta = getScheduledSetpointWeight(P_ref);
    float weightedError = P - beta * P_ref;

    _previousPressure = P;

    float epsilon = fmaxf(_epsilonCoefficient * _filteredSetpoint, 1e-4f);
    float deadband = _deadbandCoefficient * _filteredSetpoint;

    float s = _convergenceGain * weightedError;
    float sat_s = 0.0f;

    if (weightedError > 0.0f) {
        float tanv = tanhf(s / epsilon - deadband * _convergenceGain / epsilon);
        sat_s = std::max(0.0f, tanv);
    } else if (weightedError < 0.0f) {
        float tanv = tanhf(s / epsilon + deadband * _convergenceGain / epsilon);
        sat_s = std::min(0.0f, tanv);
    }

    float pressureRatio = 0.0f;
    if (P < _maxPressure) {
        pressureRatio = P / _maxPressure;
    }

    float denominator = fmaxf(1.0f - pressureRatio, 0.0001f);
    float Ki = _integralGain / denominator;

    float Qa = fmaxf(getAvailableFlow(), 1e-3f);
    float Ceq = _systemCompliance;
    float K = _commutationGain / denominator * Qa / Ceq;

    // Feedback term without integral.
    float feedbackNoIntegral = Ceq / Qa * (-_convergenceGain * weightedError - K * sat_s);

    // Dynamic FF only.
    float ffHoldPct = getPressureFeedforwardHold(P_ref);
    float ffDynRawPct = getPressureFeedforwardDynamicUnfiltered(P_ref, _filteredSetpointDerivative, error);

    // First, low-pass the dynamic FF request.
    applyLowPassFilter(&_filteredFeedforwardDynamicOutput, ffDynRawPct, _feedforwardDynamicFilterFreq, _dt);

    // Then, apply a dedicated slew limiter to the dynamic FF path itself.
    float ffDynPct = _filteredFeedforwardDynamicOutput;

    const float maxFfRisePerCycle = _feedforwardDynamicRiseRate * _dt;
    const float maxFfDropPerCycle = _feedforwardDynamicDropRate * _dt;

    if (ffDynPct > _lastFeedforwardDynamicAppliedOutput + maxFfRisePerCycle) {
        ffDynPct = _lastFeedforwardDynamicAppliedOutput + maxFfRisePerCycle;
    } else if (ffDynPct < _lastFeedforwardDynamicAppliedOutput - maxFfDropPerCycle) {
        ffDynPct = _lastFeedforwardDynamicAppliedOutput - maxFfDropPerCycle;
    }

    ffDynPct = std::clamp(ffDynPct, 0.0f, _feedforwardDynamicMaxPct);
    _lastFeedforwardDynamicAppliedOutput = ffDynPct;

    float ffTotalPct = ffHoldPct + ffDynPct;

    // Candidate integral state.
    float integralCandidate = _errorIntegral + error * _dt;

    // Candidate total output before and after 0..100 clamp.
    float feedbackPctCandidate = (feedbackNoIntegral - Ki * integralCandidate) * 100.0f;
    float totalUnclampedCandidate = feedbackPctCandidate + ffTotalPct;
    float totalClampedCandidate = std::clamp(totalUnclampedCandidate, 0.0f, 100.0f);

    bool clampHighCandidate = totalUnclampedCandidate > 100.0f;
    bool clampLowCandidate = totalUnclampedCandidate < 0.0f;

    const float maxRisePerCycle = _maxPressureOutputRiseRate * _dt;
    float maxDropPerCycle = _maxPressureOutputDropRate * _dt;

    bool limiterUpCandidate = totalClampedCandidate > (_lastAppliedPressureOutput + maxRisePerCycle);
    bool limiterDownCandidate = totalClampedCandidate < (_lastAppliedPressureOutput - maxDropPerCycle);

    // Robust anti-windup by conditional integration / clamping.
    // We stop integrating when the controller is trying to push further into a
    // real downstream limitation (hard clamp or slew limit).
    //
    // error < 0 => below target => controller wants MORE output
    // error > 0 => above target => controller wants LESS output
    bool pushingMoreOutput = error < 0.0f;
    bool pushingLessOutput = error > 0.0f;

    bool blockIntegral =
        (pushingMoreOutput && (clampHighCandidate || limiterUpCandidate)) ||
        (pushingLessOutput && (clampLowCandidate || limiterDownCandidate));

    if (!blockIntegral) {
        _errorIntegral = integralCandidate;
    }

    // Final output with accepted integral state.
    float feedbackPct = (feedbackNoIntegral - Ki * _errorIntegral) * 100.0f;
    float totalUnclamped = feedbackPct + ffTotalPct;
    float totalClamped = std::clamp(totalUnclamped, 0.0f, 100.0f);

    // Extended telemetry snapshots for CSV logging.
    _lastFeedbackOutput = feedbackPct;
    _lastFeedforwardHoldOutput = ffHoldPct;
    _lastFeedforwardDynamicOutput = ffDynPct;
    _lastFeedforwardTotalOutput = ffTotalPct;

    _lastUnclampedPressureOutput = totalUnclamped;
    _lastClampedPressureOutput = totalClamped;

    _pressureClampActiveHigh = totalUnclamped > 100.0f;
    _pressureClampActiveLow = totalUnclamped < 0.0f;

    // Backward-compatible meaning for existing telemetry:
    _lastRawPressureOutput = totalClamped;

    _pressureLimiterActiveUp = false;
    _pressureLimiterActiveDown = false;

    float appliedOutput = totalClamped;

    // Bidirectional slew limiter. This acts as a downstream actuator guardrail.
    if (totalClamped > _lastAppliedPressureOutput + maxRisePerCycle) {
        appliedOutput = _lastAppliedPressureOutput + maxRisePerCycle;
        _pressureLimiterActiveUp = true;
    } else if (totalClamped < _lastAppliedPressureOutput - maxDropPerCycle) {
        appliedOutput = _lastAppliedPressureOutput - maxDropPerCycle;
        _pressureLimiterActiveDown = true;
    }

    appliedOutput = std::clamp(appliedOutput, 0.0f, 100.0f);
    _lastAppliedPressureOutput = appliedOutput;

    // Export normalized internal pressure-branch command for telemetry.
    _pumpDutyCycle = totalClamped / 100.0f;

    return appliedOutput;
}

void PressureController::reset() {
    initSetpointFilter(_filteredPressureSensor);
    _errorIntegral = 0.0f;
    _pumpDutyCycle = 0.0f;

    // Reset all telemetry / FF / limiter state so the next shot does not inherit stale state.
    _lastFeedbackOutput = 0.0f;
    _lastFeedforwardHoldOutput = 0.0f;
    _lastFeedforwardDynamicOutput = 0.0f;
    _lastFeedforwardTotalOutput = 0.0f;
    _filteredFeedforwardDynamicOutput = 0.0f;
    _lastFeedforwardDynamicAppliedOutput = 0.0f;

    _lastUnclampedPressureOutput = 0.0f;
    _lastClampedPressureOutput = 0.0f;
    _lastRawPressureOutput = 0.0f;
    _lastAppliedPressureOutput = 0.0f;

    _pressureClampActiveHigh = false;
    _pressureClampActiveLow = false;
    _pressureLimiterActiveUp = false;
    _pressureLimiterActiveDown = false;

    _pumpFlowRate = 0.0f;
    _puckSaturationVolume = 0.0f;
    _puckState[0] = false;
    _puckState[1] = false;
    _puckState[2] = false;
    _puckCounter = 0;

    ESP_LOGI("", "RESET");
}
