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

        # Scale for higher dimensions
        scale = :math.pow(:math.sqrt(2), dims - 2)
        {scale_recursive(approx, scale), scale_recursive(details, scale)}
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
        scale = :math.pow(:math.sqrt(2), -(dims - 2))
        descaled_approx = scale_recursive(approximation, scale)
        descaled_details = scale_recursive(details, scale)

        # Reconstruct recursively
        Enum.zip(descaled_approx, descaled_details)
        |> Enum.map(fn {a, d} -> inverse(a, d, filter, dims - 1, precision) end)
    end
  end

  # Recursive scaling that handles all data types
  defp scale_recursive(data, scale) when is_number(data), do: data * scale

  defp scale_recursive(data, scale) when is_list(data) do
    Enum.map(data, &scale_recursive(&1, scale))
  end

  defp scale_recursive({a, b, c}, scale) do
    {scale_recursive(a, scale), scale_recursive(b, scale),
     scale_recursive(c, scale)}
  end

  defp scale_recursive({a, b}, scale) do
    {scale_recursive(a, scale), scale_recursive(b, scale)}
  end
end
