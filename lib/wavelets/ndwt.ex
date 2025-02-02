defmodule Wavelets.NDWT do
  @moduledoc """
  Implementation of n-dimensional Discrete Wavelet Transform without additional scaling
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

      _ ->
        transform_nd(signal, filter, dims, precision)
    end
  end

  defp transform_nd(signal, filter, dims, precision) do
    signal
    |> Enum.map(&forward(&1, filter, dims - 1, precision))
    |> Enum.unzip()
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

      _ ->
        inverse_transform_nd(approximation, details, filter, dims, precision)
    end
  end

  defp inverse_transform_nd(approximation, details, filter, dims, precision) do
    Enum.zip(approximation, details)
    |> Enum.map(fn {a, d} -> inverse(a, d, filter, dims - 1, precision) end)
  end
end
