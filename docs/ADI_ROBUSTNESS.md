# Synthetic ADI robustness diagnostics

Two independent checks complement the existing pedagogical ADI demo:

1. `adi_reduce` forms a median reference from other frames separated by more than a user-specified parallactic-angle exclusion. This reduces the most direct companion self-subtraction but does **not** remove all overlap or guarantee unbiased throughput.
2. `matched_filter_amplitude` estimates a template amplitude and formal variance using known independent pixel variances. It does not account for correlated speckles, reference-subtraction covariance, or look-elsewhere effects.

Run `python -m pytest -q tests/test_adi_robustness.py`. The tests use synthetic static speckles and exact linear models; no actual companion detections are inferred.

Authorship: repository owner Biswajit Jana. Implementation developed with AI coding assistance; scientific validation and review remain necessary.
