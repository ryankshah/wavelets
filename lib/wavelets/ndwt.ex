defmodule Wavelets.NDWT do
  @moduledoc """
  Implementation of n-dimensional Discrete Wavelet Transform
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
        # Handle higher dimensions recursively
        {approx, details} =
          signal
          |> Enum.map(&forward(&1, filter, dims - 1, precision))
          |> Enum.unzip()

        # Apply scaling for each additional dimension
        scale = :math.pow(2, (dims - 2) / 2)

        {
          deep_scale(approx, scale),
          deep_scale(details, scale)
        }
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
        # Remove dimension scaling before reconstruction
        scale = :math.pow(2, -(dims - 2) / 2)
        descaled_approx = deep_scale(approximation, scale)
        descaled_details = deep_scale(details, scale)

        # Reconstruct recursively
        Enum.zip(descaled_approx, descaled_details)
        |> Enum.map(fn {a, d} -> inverse(a, d, filter, dims - 1, precision) end)
    end
  end

  defp deep_scale(data, scale) when is_number(data), do: data * scale

  defp deep_scale(data, scale) when is_list(data) do
    Enum.map(data, &deep_scale(&1, scale))
  end
end
