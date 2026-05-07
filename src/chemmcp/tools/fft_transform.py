import logging
import math
import cmath
from typing import Dict, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


def _dft(signal: List[complex]) -> List[complex]:
    """Pure Python DFT (slow but correct for fallback)."""
    N = len(signal)
    result = []
    for k in range(N):
        val = complex(0, 0)
        for n in range(N):
            angle = 2 * math.pi * k * n / N
            val += signal[n] * complex(math.cos(angle), -math.sin(angle))
        result.append(val / N)  # normalize
    return result


def _apply_window(signal: List[float], window: str) -> List[float]:
    """Apply window function to signal."""
    N = len(signal)
    if window == "none" or window is None:
        return list(signal)

    w = [1.0] * N
    if window == "hanning":
        for i in range(N):
            w[i] = 0.5 * (1 - math.cos(2 * math.pi * i / (N - 1)))
    elif window == "hamming":
        for i in range(N):
            w[i] = 0.54 - 0.46 * math.cos(2 * math.pi * i / (N - 1))
    elif window == "blackman":
        for i in range(N):
            w[i] = 0.42 - 0.5 * math.cos(2 * math.pi * i / (N - 1)) + \
                   0.08 * math.cos(4 * math.pi * i / (N - 1))
    else:
        raise ChemMCPError(f"Unknown window type: {window}. Use: none/hanning/hamming/blackman")

    return [signal[i] * w[i] for i in range(N)]


@ChemMCPManager.register_tool
class FftTransform(BaseTool):
    """
    快速傅里叶变换工具 —— 光谱处理、NMR分析。
    对时域信号进行FFT变换，得到频谱信息。
    """
    __version__ = "0.1.0"
    name = "FftTransform"
    func_name = "fft_transform"
    description = (
        "Perform Fast Fourier Transform (FFT) on time-domain signal data. "
        "Essential for spectroscopy processing, NMR analysis, and frequency-domain analysis."
    )
    implementation_description = (
        "Uses numpy.fft when available for fast computation, "
        "otherwise falls back to pure Python DFT. Supports Hanning/Hamming/Blackman windows."
    )
    oss_dependencies = [
        ("numpy", "https://numpy.org/", "BSD"),
    ]
    services_and_software = []
    categories = ["General"]
    tags = ["FFT", "Spectroscopy", "Signal Processing", "NMR", "Frequency Analysis"]
    required_envs = []

    code_input_sig = [
        ("signal_data", "list", "N/A", "Time-domain signal data (list of floats or complex numbers)."),
        ("sampling_rate", "float", "N/A", "Sampling rate in Hz."),
        ("window", "str", "none", "Window function: 'none', 'hanning', 'hamming', or 'blackman'."),
        ("zero_padding", "int", "0", "Zero-padding length to improve frequency resolution (default: 0)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A",
         "Format: 'v1,v2,v3,...; sampling_rate; [window]; [zero_padding]' "
         "Example: 'sin data values; 1000; hanning; 0'"),
    ]

    output_sig = [
        ("frequencies", "list", "Frequency axis in Hz (positive frequencies only)."),
        ("amplitudes", "list", "Amplitude spectrum (magnitude at each frequency)."),
        ("phases", "list", "Phase spectrum in radians at each frequency."),
        ("dominant_frequency", "float", "Frequency with maximum amplitude."),
        ("dominant_amplitude", "float", "Maximum amplitude value."),
    ]

    examples = [
        {
            "code_input": {
                "signal_data": [math.sin(2 * math.pi * 50 * t / 1000) + 0.5 * math.sin(2 * math.pi * 120 * t / 1000)
                                 for t in range(256)],
                "sampling_rate": 1000.0,
                "window": "none",
                "zero_padding": 0,
            },
            "text_input": {
                "input_str": "generated_sine; 1000; none; 0",
            },
            "output": {
                "dominant_frequency": 50.0,
                "dominant_amplitude": 0.5,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._use_numpy = False
        try:
            import numpy as np
            self._numpy = np
            self._use_numpy = True
        except ImportError:
            logger.info("numpy not available, using DFT fallback")

    def _run_base(
        self,
        signal_data: List,
        sampling_rate: float,
        window: str = "none",
        zero_padding: int = 0,
    ) -> Dict:
        if len(signal_data) < 2:
            raise ChemMCPError("Signal must have at least 2 data points.")

        # Convert to float list
        try:
            signal = [complex(s) if isinstance(s, complex) else float(s) for s in signal_data]
        except (ValueError, TypeError) as e:
            raise ChemMCPError(f"Invalid signal data: {e}")

        # Apply window
        real_signal = [s.real for s in signal]
        windowed = _apply_window(real_signal, window.lower())
        signal = [complex(w) for w in windowed]

        # Zero padding
        if zero_padding > 0:
            signal = signal + [complex(0)] * zero_padding

        N = len(signal)

        # FFT
        if self._use_numpy:
            fft_result = self._numpy.fft.fft(signal) / N
        else:
            fft_result = _dft(signal)

        # Positive frequencies only
        n_pos = N // 2
        freqs = [(k * sampling_rate / N) for k in range(n_pos)]
        amplitudes = [abs(fft_result[k]) for k in range(n_pos)]
        phases = [cmath.phase(fft_result[k]) for k in range(n_pos)]

        # Find dominant frequency
        max_idx = 0
        max_amp = amplitudes[0] if amplitudes else 0
        for i in range(1, len(amplitudes)):
            if amplitudes[i] > max_amp:
                max_amp = amplitudes[i]
                max_idx = i

        result = {
            "frequencies": [round(f, 6) for f in freqs],
            "amplitudes": [round(a, 6) for a in amplitudes],
            "phases": [round(p, 6) for p in phases],
            "dominant_frequency": round(freqs[max_idx], 4) if freqs else 0.0,
            "dominant_amplitude": round(max_amp, 6),
        }
        logger.info(f"FFT done: N={N}, fs={sampling_rate}Hz, dominant_freq={result['dominant_frequency']}Hz")
        return result

    def _run_text(self, input_str: str) -> Dict:
        try:
            parts = [p.strip() for p in input_str.split(";")]
            if len(parts) < 2:
                raise ValueError("Need at least 2 parts: signal_data; sampling_rate")

            signal_data = [float(v) for v in parts[0].split(",")]
            sampling_rate = float(parts[1])
            window = parts[2].strip().lower() if len(parts) > 2 else "none"
            zero_padding = int(parts[3]) if len(parts) > 3 else 0

            return self._run_base(signal_data, sampling_rate, window, zero_padding)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
