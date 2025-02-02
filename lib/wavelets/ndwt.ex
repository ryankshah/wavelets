defmodule Wavelets.NDWT do
  @moduledoc """
  Implementation of n-dimensional Discrete Wavelet Transform with consistent scaling
  """

  alias Wavelets.{DWT, Filter, IDWT}

  @doc """
  Performs n-dimensional forward discrete wavelet transform
  """
  def forward(signal, filter, dims, precision \\ :double) when dims > 0 do
    case dims do
      1 when is_list(signal) ->
        DWT.forward_1d(signal, filter, precision)

      2 when is_list(hd(signal)) ->
        DWT.forward_2d(signal, filter, precision)

      _ when dims > 2 ->
        # Apply transform recursively for higher dimensions
        signal
        |> Enum.map(&forward(&1, filter, dims - 1, precision))
        |> Enum.unzip()
    end
  end

  @doc """
  Performs n-dimensional inverse discrete wavelet transform
  """
  def inverse(approximation, details, filter, dims, precision \\ :double)
      when dims > 0 do
    case dims do
      1 ->
        IDWT.inverse_1d(approximation, details, filter, precision)

      2 ->
        IDWT.inverse_2d(approximation, details, filter, precision)

      _ when dims > 2 ->
        # Apply inverse transform recursively
        Enum.zip(approximation, details)
        |> Enum.map(fn {a, d} -> inverse(a, d, filter, dims - 1, precision) end)
    end
  end
end
