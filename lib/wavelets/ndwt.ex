defmodule Wavelets.NDWT do
  @moduledoc """
  Implementation of n-dimensional Discrete Wavelet Transform with simple scaling
  """

  alias Wavelets.{DWT, Filter, IDWT}

  @doc """
  Performs n-dimensional forward discrete wavelet transform
  """
  def forward(signal, filter, dims, precision \\ :double) when dims > 0 do
    # Scale input by 2.0 for each dimension
    scale = :math.pow(2.0, dims)
    scaled_signal = scale_nd(signal, scale, dims)

    case dims do
      1 when is_list(signal) ->
        DWT.forward_1d(scaled_signal, filter, precision)

      2 when is_list(hd(signal)) ->
        DWT.forward_2d(scaled_signal, filter, precision)

      _ ->
        transform_nd(scaled_signal, filter, dims, precision)
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
    result =
      case dims do
        1 ->
          IDWT.inverse_1d(approximation, details, filter, precision)

        2 ->
          IDWT.inverse_2d(approximation, details, filter, precision)

        _ ->
          inverse_transform_nd(approximation, details, filter, dims, precision)
      end

    # Unscale result by 2.0 for each dimension
    scale = 1 / :math.pow(2.0, dims)
    scale_nd(result, scale, dims)
  end

  defp inverse_transform_nd(approximation, details, filter, dims, precision) do
    Enum.zip(approximation, details)
    |> Enum.map(fn {a, d} -> inverse(a, d, filter, dims - 1, precision) end)
  end

  defp scale_nd(data, scale, dims) when dims > 0 do
    case data do
      list when is_list(list) ->
        if dims == 1 do
          Enum.map(list, &(&1 * scale))
        else
          Enum.map(list, &scale_nd(&1, scale, dims - 1))
        end

      value ->
        value * scale
    end
  end
end
