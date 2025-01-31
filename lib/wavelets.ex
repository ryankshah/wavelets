defmodule Wavelets do
  @moduledoc """
  A comprehensive wavelet transform library for Elixir.
  Provides implementations of various wavelet transforms including:
  - Discrete Wavelet Transform (DWT)
  - Inverse Discrete Wavelet Transform (IDWT)
  - Stationary Wavelet Transform
  - Wavelet Packet Transform
  - Continuous Wavelet Transform

  Supports both real and complex calculations in single and double precision.
  """

  @type precision :: :single | :double
  @type complex :: {float, float}
end
