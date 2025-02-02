defmodule Wavelets.NDWT do
  @moduledoc """
  Implementation of n-dimensional Discrete Wavelet Transform with dimension-appropriate scaling
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
        # Adjust scaling for higher dimensions
        {approx, details} = transform_nd(signal, filter, dims, precision)
        scale = :math.pow(:math.sqrt(2), dims - 2)

        {
          scale_nd(approx, scale, dims),
          scale_nd(details, scale, dims)
        }
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
        # Adjust scaling for inverse transform
        scale = :math.pow(:math.sqrt(2), 2 - dims)
        scaled_approx = scale_nd(approximation, scale, dims)
        scaled_details = scale_nd(details, scale, dims)

        inverse_transform_nd(
          scaled_approx,
          scaled_details,
          filter,
          dims,
          precision
        )
    end
  end

  defp inverse_transform_nd(approximation, details, filter, dims, precision) do
    Enum.zip(approximation, details)
    |> Enum.map(fn {a, d} -> inverse(a, d, filter, dims - 1, precision) end)
  end

  defp scale_nd(data, scale, dims) when dims > 2 do
    case data do
      list when is_list(list) ->
        Enum.map(list, &scale_nd(&1, scale, dims - 1))

      value ->
        value * scale
    end
  end

  defp scale_nd(data, _scale, _dims), do: data
end
