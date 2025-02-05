defmodule Wavelets.NDWT do
  @moduledoc """
  Implementation of n-dimensional Discrete Wavelet Transform with proper normalization
  """

  alias Wavelets.{DWT, Filter, IDWT}

  @doc """
  Performs n-dimensional forward discrete wavelet transform
  """
  def forward(signal, filter, dims, precision \\ :double) when dims > 0 do
    # Simple recursive application of 1D transform
    case dims do
      1 when is_list(signal) ->
        DWT.forward_1d(signal, filter, precision)

      2 when is_list(hd(signal)) ->
        DWT.forward_2d(signal, filter, precision)

      _ when dims > 2 ->
        signal
        |> Enum.map(&forward(&1, filter, dims - 1, precision))
        |> Enum.unzip()
    end
  end

  defp scale_signal(data, scale) when is_number(data), do: data * scale

  defp scale_signal(data, scale) when is_list(data) do
    Enum.map(data, &scale_signal(&1, scale))
  end

  defp scale_signal({a, b}, scale) do
    {scale_signal(a, scale), scale_signal(b, scale)}
  end

  defp scale_signal({a, b, c}, scale) do
    {scale_signal(a, scale), scale_signal(b, scale), scale_signal(c, scale)}
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
        # Inverse transform recursively
        Enum.zip(approximation, details)
        |> Enum.map(fn {a, d} ->
          inverse(a, d, filter, dims - 1, precision)
        end)
    end
  end

  # Recursive scaling functions handling all data types
  defp scale_recursive(data, scale) when is_number(data), do: data * scale

  defp scale_recursive(data, scale) when is_list(data) do
    Enum.map(data, &scale_recursive(&1, scale))
  end

  defp scale_recursive({a, b, c}, scale) do
    {scale_recursive(a, scale), scale_recursive(b, scale),
     scale_recursive(c, scale)}
  end
end
